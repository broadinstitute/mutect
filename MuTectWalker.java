package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import java.io.*;
import java.util.*;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.*;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BetaDistribution;
import org.apache.commons.math.distribution.BetaDistributionImpl;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotype;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

@PartitionBy(PartitionType.LOCUS)
@BAQMode()
@Reference(window=@Window(start=-1* MuTectWalker.REFERENCE_HALF_WINDOW_LENGTH,stop= MuTectWalker.REFERENCE_HALF_WINDOW_LENGTH))
@By(DataSource.REFERENCE)
public class MuTectWalker extends LocusWalker<Integer, Integer> implements TreeReducible<Integer> {
    public static final int REFERENCE_HALF_WINDOW_LENGTH = 150;
    public static final String BAM_TAG_TUMOR = "tumor";
    public static final String BAM_TAG_NORMAL = "normal";

    // DO NOT CHANGE THIS LINE!  It's the SVN revision number of the caller, which updates automatically!
    private static final String VERSION = "$Rev$";

    @ArgumentCollection private MuTectArgumentCollection MTAC = new MuTectArgumentCollection();

    /***************************************/
    // Call-stats output
    /***************************************/
    @Output(doc="Call-stats output")
    PrintStream out;

    /***************************************/
    // Reference Metadata inputs
    /***************************************/
    @Input(fullName="dbsnp", shortName = "dbsnp", doc="VCF file of DBSNP information", required=false)
    public RodBinding<VariantContext> dbsnpRod;

    @Input(fullName="cosmic", shortName = "cosmic", doc="VCF file of COSMIC sites", required=false)
    public RodBinding<VariantContext> cosmicRod;

    @Hidden
    @Input(fullName="normal_panel", shortName = "normal_panel", doc="VCF file of sites observed in normal", required=false)
    public RodBinding<VariantContext> normalPanelRod;

    /***************************************/
    // coverage outputs
    /***************************************/
    @Output(fullName="coverage_file", shortName="cov", doc="write out coverage in WIGGLE format to this file", required=false)
    public PrintStream COVERAGE_FILE = null;

    @Output(fullName="coverage_20_q20_file", shortName="cov_q20", doc="write out 20x of Q20 coverage in WIGGLE format to this file", required=false)
    public PrintStream COVERAGE_20_Q20_FILE = null;

    @Output(fullName="power_file", shortName="pow", doc="write out power in WIGGLE format to this file", required=false)
    public PrintStream POWER_FILE = null;

    @Output(fullName="tumor_depth_file", shortName="tdf", doc="write out tumor read depth in WIGGLE format to this file", required=false)
    public PrintStream TUMOR_DEPTH_FILE = null;

    @Output(fullName="normal_depth_file", shortName="ndf", doc="write out normal read depth in WIGGLE format to this file", required=false)
    public PrintStream NORMAL_DEPTH_FILE = null;

    public boolean NO_BAQ = true;
    public int MIN_QSUM_QSCORE = 13;
    public boolean USE_MAPQ0_IN_NORMAL_QSCORE = true;

    private static final FisherExact fisher = new FisherExact(5000);

    private boolean hasTumorBam = false;
    private boolean hasNormalBam = false;

    private double contaminantAlternateFraction;

    private double powerConstantEps;
    private TumorPowerCalculator tumorPowerCalculator;
    private NormalPowerCalculator normalNovelSitePowerCalculator;
    private NormalPowerCalculator normalDbSNPSitePowerCalculator;
    private TumorPowerCalculator normalArtifactPowerCalculator;
    private TumorPowerCalculator strandArtifactPowerCalculator;

    private static class ThreadLocalRankSumTest extends ThreadLocal<RankSumTest> {
        public RankSumTest initialValue() {
            return new RankSumTest();
        }
    }

    private static class PileupComparatorByAltRefQual implements Comparator<PileupElement> {
        private byte ref;
        private byte alt;

        private PileupComparatorByAltRefQual(byte ref, byte alt) {
            this.ref = ref;
            this.alt = alt;
        }

        public int compare(PileupElement o1, PileupElement o2) {
            // if the bases are the same, the higher quality score comes first
            if (o1.getBase() ==  o2.getBase()) {
                if (o1.getQual() == o2.getQual()) { return 0; }
                return (o1.getQual() > o2.getQual())?-1:1;
            } else {
                return (o1.getBase() == alt)?-1:1;
            }

        }
    }


    private static ThreadLocalRankSumTest __rankSumTest = new ThreadLocalRankSumTest();

    public static RankSumTest getRankSumTest() {
        return __rankSumTest.get();
    }

    private CoverageWiggleFileWriter stdCovWriter;
    private CoverageWiggleFileWriter q20CovWriter;
    private CoverageWiggleFileWriter powerWriter;
    private CoverageWiggleFileWriter tumorDepthWriter;
    private CoverageWiggleFileWriter normalDepthWriter;


    private Set<SAMReaderID> tumorSAMReaderIDs = new HashSet<SAMReaderID>();

    private CallStatsGenerator callStatsGenerator;

    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    @Override
	public void initialize() {
        if (MTAC.NOOP) { return; }

        refReader = this.getToolkit().getReferenceDataSource().getReference();
        callStatsGenerator = new CallStatsGenerator(MTAC.ENABLE_EXTENDED_OUTPUT);

        // check that we have at least one tumor bam
        for(SAMReaderID id : getToolkit().getReadsDataSource().getReaderIDs()) {
            if (id.getTags().getPositionalTags().size() == 0) {
                throw new RuntimeException("BAMs must be tagged as either 'tumor' or 'normal'");
            }

            for(String tag : id.getTags().getPositionalTags()) {
                if (BAM_TAG_TUMOR.equalsIgnoreCase(tag)) {
                    hasTumorBam = true;
                    tumorSAMReaderIDs.add(id);

                    // fill in the sample name if necessary
                    if (MTAC.TUMOR_SAMPLE_NAME == null) {
                        try {
                            if (getToolkit().getReadsDataSource().getHeader(id).getReadGroups().size() == 0) {
                                throw new RuntimeException("No Read Groups found for Tumor BAM -- Read Groups are Required, or supply tumor_sample_name!");
                            }
                            MTAC.TUMOR_SAMPLE_NAME = getToolkit().getReadsDataSource().getHeader(id).getReadGroups().get(0).getSample();
                        } catch (NullPointerException npe) {
                            MTAC.TUMOR_SAMPLE_NAME = "tumor";
                        }
                    }
                } else if (BAM_TAG_NORMAL.equalsIgnoreCase(tag)) {
                    hasNormalBam = true;

                    // fill in the sample name if necessary
                    if (MTAC.NORMAL_SAMPLE_NAME == null) {
                        try {
                            if (getToolkit().getReadsDataSource().getHeader(id).getReadGroups().size() == 0) {
                                throw new RuntimeException("No Read Groups found for Normal BAM -- Read Groups are Required, or supply normal_sample_name!");
                            }

                            MTAC.NORMAL_SAMPLE_NAME = getToolkit().getReadsDataSource().getHeader(id).getReadGroups().get(0).getSample();
                        } catch (NullPointerException npe) {
                            MTAC.NORMAL_SAMPLE_NAME = "normal";
                        }
                    }
                } else {
                    throw new RuntimeException("Unknown BAM tag '" + tag + "' must be either 'tumor' or 'normal'");
                }                
            }
        }

        if (!hasTumorBam) {
            throw new RuntimeException("At least one BAM tagged as 'tumor' required");
        }

        if (!hasNormalBam) {
            MTAC.NORMAL_LOD_THRESHOLD = -1 * Float.MAX_VALUE;
            MTAC.NORMAL_DBSNP_LOD_THRESHOLD = -1 * Float.MAX_VALUE;
            MTAC.NORMAL_ARTIFACT_LOD_THRESHOLD = Float.MAX_VALUE;
        }

        this.contaminantAlternateFraction = Math.max(MTAC.MINIMUM_MUTATION_CELL_FRACTION, MTAC.FRACTION_CONTAMINATION);

        // coverage related initialization
        this.powerConstantEps = Math.pow(10, -1 * (MTAC.POWER_CONSTANT_QSCORE/10));

        this.tumorPowerCalculator = new TumorPowerCalculator(this.powerConstantEps, MTAC.TUMOR_LOD_THRESHOLD, this.contaminantAlternateFraction);
        this.normalNovelSitePowerCalculator = new NormalPowerCalculator(this.powerConstantEps, MTAC.NORMAL_LOD_THRESHOLD);
        this.normalDbSNPSitePowerCalculator = new NormalPowerCalculator(this.powerConstantEps, MTAC.NORMAL_DBSNP_LOD_THRESHOLD);
        this.normalArtifactPowerCalculator = new TumorPowerCalculator(this.powerConstantEps, MTAC.NORMAL_ARTIFACT_LOD_THRESHOLD, 0.0f);
        this.strandArtifactPowerCalculator = new TumorPowerCalculator(this.powerConstantEps, MTAC.STRAND_ARTIFACT_LOD_THRESHOLD, 0.0f);

        stdCovWriter = new CoverageWiggleFileWriter(COVERAGE_FILE);
        q20CovWriter = new CoverageWiggleFileWriter(COVERAGE_20_Q20_FILE);
        powerWriter = new CoverageWiggleFileWriter(POWER_FILE);
        tumorDepthWriter = new CoverageWiggleFileWriter(TUMOR_DEPTH_FILE);
        normalDepthWriter = new CoverageWiggleFileWriter(NORMAL_DEPTH_FILE);

        // to force output, all we have to do is lower the initial tumor lod threshold to -infinity
        if (MTAC.FORCE_OUTPUT) {
            MTAC.INITIAL_TUMOR_LOD_THRESHOLD = -Float.MAX_VALUE;
        }

        // write out the call stats header
        out.println("## muTector v1.0." + VERSION.split(" ")[1]);
        out.println(callStatsGenerator.generateHeader());

        lastTime = System.currentTimeMillis();
    }

    public static int MAX_INSERT_SIZE = 10000;
    private int totalReadsProcessed = 0;
    private int binReadsProcessed = 0;
    private long lastTime;
    private int candidatesInspected = 0;

    @Override
	public Integer map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext rawContext) {
        if (MTAC.NOOP) return 0;
        
        TreeMap<Double, String> messageByTumorLod = new TreeMap<Double, String>();

        ReadBackedPileup pileup = rawContext.getBasePileup();
        int numberOfReads = pileup.depthOfCoverage();
        binReadsProcessed += numberOfReads;

        if (binReadsProcessed >= 1000000) {
            long time = System.currentTimeMillis();
            long elapsedTime = time - lastTime;
            lastTime = time;

            totalReadsProcessed += binReadsProcessed;
            binReadsProcessed = 0;

            logger.info(String.format("[MUTECTOR] Processed %d reads in %d ms", totalReadsProcessed, elapsedTime));
        }

        // an optimization to speed things up when there is no coverage
        if ( !MTAC.FORCE_OUTPUT && numberOfReads == 0) { return -1; }

        final char upRef = Character.toUpperCase(ref.getBaseAsChar());
        String sequenceContext = createSequenceContext(ref, 3);

        try {


            // only process bases where the reference is [ACGT], because the FASTA for HG18 has N,M and R!
            if (upRef != 'A' && upRef != 'C' && upRef != 'G' && upRef != 'T') {
                return -1;
            }

            ArrayList<PileupElement> normalPileupElements = new ArrayList<PileupElement>();
            ArrayList<PileupElement> tumorPileupElements = new ArrayList<PileupElement>();

            int totalPairs = 0;
            int improperPairs = 0;
            for (PileupElement p : pileup ) {
                final GATKSAMRecord read = p.getRead();
                final byte base = p.getBase();

                if (base == ((byte)'N') || base == ((byte)'n')) {
                    continue;
                }

                // count # of reads that are part of a mapped pair where the 'proper pair'
                // flag is false or the insert size is > 10kb
                if (read.getReadPairedFlag() && !read.getMateUnmappedFlag()) {
                    totalPairs++;
                    if (!read.getProperPairFlag() || read.getInferredInsertSize() >= MAX_INSERT_SIZE) {
                        improperPairs++;
                    }
                }


                // Add the read to the appropriate pile of reads
                if (isTumorRead(p.getRead())) {
                    tumorPileupElements.add(p);
                } else {
                    normalPileupElements.add(p);
                }
            }

            ReadBackedPileup normalPileup =
                    new ReadBackedPileupImpl(rawContext.getLocation(), normalPileupElements);

            ReadBackedPileup tumorPileup =
                    new ReadBackedPileupImpl(rawContext.getLocation(), tumorPileupElements);

            if (MTAC.BAM_TUMOR_SAMPLE_NAME != null) {
                tumorPileup = tumorPileup.getPileupForSample(MTAC.BAM_TUMOR_SAMPLE_NAME);
            }

            final LocusReadPile tumorReadPile = new LocusReadPile(tumorPileup, upRef, MTAC.MIN_QSCORE, MIN_QSUM_QSCORE, false, MTAC.ARTIFACT_DETECTION_MODE);
            final LocusReadPile normalReadPile = new LocusReadPile(normalPileup, upRef, MTAC.MIN_QSCORE, 0, this.USE_MAPQ0_IN_NORMAL_QSCORE, true);


            Collection<VariantContext> panelOfNormalsVC = tracker.getValues(normalPanelRod, rawContext.getLocation());
            Collection<VariantContext> cosmicVC = tracker.getValues(cosmicRod, rawContext.getLocation());
            Collection<VariantContext> dbsnpVC = tracker.getValues(dbsnpRod, rawContext.getLocation());

            // remove the effect of cosmic from dbSNP
            boolean germlineAtRisk = (!dbsnpVC.isEmpty() && cosmicVC.isEmpty());

            // compute coverage flags
            int tumorCoveredDepthThreshold = 14;
            int normalCoveredDepthThreshold = (germlineAtRisk)?19:8;
            if (!hasNormalBam) {
                normalCoveredDepthThreshold = 0;
            }

            int tumorBaseCount = tumorReadPile.finalPileupReads.size();
            int normalBaseCount = normalReadPile.finalPileupReads.size();
            boolean isTumorCovered = tumorBaseCount >= tumorCoveredDepthThreshold;
            boolean isNormalCovered = normalBaseCount >= normalCoveredDepthThreshold;
            boolean isBaseCovered = isTumorCovered && isNormalCovered;
            if (!hasNormalBam) {
                isBaseCovered = isTumorCovered;
            }

            stdCovWriter.writeCoverage(rawContext, isBaseCovered);
            int tumorQ20BaseCount = tumorReadPile.getFilteredBaseCount(20);
            int normalQ20BaseCount = normalReadPile.getFilteredBaseCount(20);
            q20CovWriter.writeCoverage(rawContext, tumorQ20BaseCount  >= 20 &&  normalQ20BaseCount >= 20);
            tumorDepthWriter.writeCoverage(rawContext, tumorBaseCount);
            normalDepthWriter.writeCoverage(rawContext, normalBaseCount);

            // calculate power
            // TODO: use ASCN information
            double tumorPower;
            double normalPower;
            double normalPowerWithSNPPrior;
            double normalPowerNoSNPPrior;
            double combinedPower;
            if (MTAC.ABSOLUTE_COPY_NUMBER_DATA == null) {
                tumorPower = tumorPowerCalculator.cachingPowerCalculation(tumorBaseCount, MTAC.POWER_CONSTANT_AF);

                normalPowerNoSNPPrior = normalNovelSitePowerCalculator.cachingPowerCalculation(normalBaseCount);
                normalPowerWithSNPPrior = normalDbSNPSitePowerCalculator.cachingPowerCalculation(normalBaseCount);
                
                normalPower = (germlineAtRisk)?normalPowerWithSNPPrior:normalPowerNoSNPPrior;
                
                combinedPower = tumorPower*normalPower;
                if (!hasNormalBam) {
                    combinedPower = tumorPower;
                }

                powerWriter.writeCoverage(rawContext, combinedPower);
            } else {
                throw new RuntimeException("ASCN Power Calculations Not Yet Implemented!");
            }

            int mapQ0Reads =
                    tumorReadPile.qualityScoreFilteredPileup.getNumberOfMappingQualityZeroReads() +
                    normalReadPile.qualityScoreFilteredPileup.getNumberOfMappingQualityZeroReads();


            // Test each of the possible alternate alleles
            for (final char altAllele : new char[]{'A','C','G','T'}) {
                if (altAllele == upRef) { continue; }
                if (!MTAC.FORCE_OUTPUT && tumorReadPile.qualitySums.getCounts(altAllele) == 0) { continue; }

                CandidateMutation candidate = new CandidateMutation(rawContext.getLocation(), upRef);
                candidate.setSequenceContext(sequenceContext);
                candidate.setTumorSampleName(MTAC.TUMOR_SAMPLE_NAME);
                candidate.setNormalSampleName(MTAC.NORMAL_SAMPLE_NAME);
                candidate.setCovered(isBaseCovered);
                candidate.setPower(combinedPower);
                candidate.setTumorPower(tumorPower);
                candidate.setNormalPower(normalPower);
                candidate.setNormalPowerWithSNPPrior(normalPowerWithSNPPrior);
                candidate.setNormalPowerNoSNPPrior(normalPowerNoSNPPrior);
                candidate.setTumorQ20Count(tumorQ20BaseCount);
                candidate.setNormalQ20Count(normalQ20BaseCount);
                candidate.setInitialTumorNonRefQualitySum(tumorReadPile.qualitySums.getOtherQualities(upRef));
                candidate.setAltAllele(altAllele);
                candidate.setTotalPairs(totalPairs);
                candidate.setImproperPairs(improperPairs);
                candidate.setMapQ0Reads(mapQ0Reads);
                candidate.setContaminationFraction(MTAC.FRACTION_CONTAMINATION);
                candidate.setPanelOfNormalsVC(panelOfNormalsVC.isEmpty()?null:panelOfNormalsVC.iterator().next());
                candidate.setCosmicSite(!cosmicVC.isEmpty());
                candidate.setDbsnpSite(!dbsnpVC.isEmpty());
                candidate.setTumorF(tumorReadPile.estimateAlleleFraction(upRef, altAllele));

                if (!MTAC.FORCE_OUTPUT && candidate.getTumorF() < MTAC.TUMOR_F_PRETEST) {
                    continue;
                }

                if (++candidatesInspected % 1000 == 0) {
                    logger.info(String.format("[MUTECTOR] Inspected %d potential candidates", candidatesInspected));
                }
                
                candidate.setInitialTumorAltCounts(tumorReadPile.qualitySums.getCounts(altAllele));
                candidate.setInitialTumorRefCounts(tumorReadPile.qualitySums.getCounts(upRef));
                candidate.setInitialTumorAltQualitySum(tumorReadPile.qualitySums.getQualitySum(altAllele));
                candidate.setInitialTumorRefQualitySum(tumorReadPile.qualitySums.getQualitySum(upRef));

                // TODO: why extract the counts twice?  once above and once in this method...
                double tumorFLB = 0;
                int refCount = candidate.getInitialTumorRefCounts();
                int altCount = candidate.getInitialTumorAltCounts();
                int depth = refCount + altCount;
                if ( altCount > 0) {
                    // implemented as shown http://www.sigmazone.com/binomial_confidence_interval.htm
                    BetaDistribution dist = new BetaDistributionImpl(depth - altCount + 1, altCount);
                    tumorFLB = 1 - dist.inverseCumulativeProbability(1 - (1-0.95)/2);
                }
                candidate.setTumorFLowerBound(tumorFLB);


                double tumorLod = tumorReadPile.calculateAltVsRefLOD((byte)altAllele, candidate.getTumorF(), 0);
                candidate.setTumorLodFStar(tumorLod);

                candidate.setInitialTumorReadDepth(tumorReadPile.finalPileupReads.size());
                candidate.setTumorInsertionCount(tumorReadPile.getInsertionsCount());
                candidate.setTumorDeletionCount(tumorReadPile.getDeletionsCount());

                if (candidate.getTumorLodFStar() < MTAC.INITIAL_TUMOR_LOD_THRESHOLD ) {
                    continue;
                }

                // calculate lod of contaminant
                double contaminantF = Math.min(contaminantAlternateFraction, candidate.getTumorF());
                VariableAllelicRatioGenotypeLikelihoods contaminantLikelihoods
                        = new VariableAllelicRatioGenotypeLikelihoods(upRef, contaminantF);

                List<PileupElement> peList = new ArrayList<PileupElement>(tumorReadPile.finalPileup.depthOfCoverage());
                for(PileupElement pe : tumorReadPile.finalPileup) {
                    peList.add(pe);
                }

                Collections.sort(peList, new PileupComparatorByAltRefQual((byte)upRef,(byte)altAllele));
                int readsToKeep = (int) (peList.size() * contaminantAlternateFraction);

                for(PileupElement pe : peList) {
                    byte base = pe.getBase();
                    if (pe.getBase() == altAllele) {
                        // if we've retained all we need, then turn the remainder of alts to ref
                        if (readsToKeep == 0) {
                            base = (byte) upRef;
                        } else {
                            readsToKeep--;
                        }
                    }

                    // TODO: when fragment-based calling is fully enabled in VARGL, this needs to change!
                    // ignore bad bases, cap at mapping quality, and effectively have no threshold on quality (80) since this was done in our pileup
                    // TODO: remove this local copy of qualToUse once we can just use the real thing!
                    // TODO: move to this                           contaminantLikelihoods.add(base, qualToUse(pe, true, true, 0));
                    contaminantLikelihoods.add(base, qualToUse(pe, false, false, 0));
                }
                double[] refHetHom = LocusReadPile.extractRefHetHom(contaminantLikelihoods, upRef, altAllele);
                double contaminantLod = refHetHom[1] - refHetHom[0];
                candidate.setContaminantLod(contaminantLod);


                // (ii) the quality score sum for the mutant base in the normal must be < 50 and the
                //      LOD score for ref:ref vs mutant:ref + mutant:mutant must be at least 2.3.
                final QualitySums normQs = normalReadPile.qualitySums;

                
                VariableAllelicRatioGenotypeLikelihoods normalGl = normalReadPile.calculateLikelihoods(normalReadPile.qualityScoreFilteredPileup); // use MAPQ0 reads
                candidate.setInitialNormalBestGenotype(normalReadPile.getBestGenotype(normalGl));
                candidate.setInitialNormalLod(LocusReadPile.getRefVsAlt(normalGl, upRef, altAllele));

// /                candidate.setInitialNormalLod(normalReadPile.calculateRefVsAltLOD(normalReadPile.qualityScoreFilteredPileup, (byte)altAllele, 0.5, 0.0));

                double normalF = Math.max(normalReadPile.estimateAlleleFraction(normalReadPile.qualityScoreFilteredPileup, upRef, altAllele), MTAC.MINIMUM_NORMAL_ALLELE_FRACTION);
                candidate.setNormalF(normalF);
                candidate.setNormalLodFStar(normalReadPile.calculateRefVsAltLOD(normalReadPile.qualityScoreFilteredPileup, (byte)altAllele, normalF, 0.0));


                // calculate power to have detected this artifact in the normal
                candidate.setNormalArtifactPowerTF(this.normalArtifactPowerCalculator.cachingPowerCalculation(normalBaseCount, candidate.getTumorF()));
                candidate.setNormalArtifactPowerLowTF(this.normalArtifactPowerCalculator.cachingPowerCalculation(normalBaseCount, candidate.getTumorFLowerBound()));
                candidate.setNormalArtifactPowerNF(this.normalArtifactPowerCalculator.cachingPowerCalculation(normalBaseCount, candidate.getNormalF()));


                // compare the local and global error models
                RecalibratedLocalQualityScores lqs = new RecalibratedLocalQualityScores((byte) upRef, normalReadPile.finalPileup);
                double lodOriginalQualities = LocusReadPile.calculateLogLikelihood(normalReadPile.finalPileup, (byte) upRef, (byte) altAllele, 0.0);
                double lodLqs = LocusReadPile.calculateLogLikelihood(normalReadPile.finalPileup, (byte) upRef, (byte) altAllele, 0.0, lqs);
                double qLod = lodLqs - lodOriginalQualities;
                candidate.setNormalGlobalQualityReferenceLL(lodOriginalQualities);
                candidate.setNormalLocalQualityReferenceLL(lodLqs);
                candidate.setNormalQualityModelLod(qLod);

                VariableAllelicRatioGenotypeLikelihoods normalArtifactGlTF = normalReadPile.calculateLikelihoods(candidate.getTumorF(), normalReadPile.qualityScoreFilteredPileup);
                candidate.setNormalArtifactLodTF(normalReadPile.getAltVsRef(normalArtifactGlTF, upRef, altAllele));

                VariableAllelicRatioGenotypeLikelihoods normalArtifactGlLowTF = normalReadPile.calculateLikelihoods(candidate.getTumorFLowerBound(), normalReadPile.qualityScoreFilteredPileup);
                candidate.setNormalArtifactLodLowTF(normalReadPile.getAltVsRef(normalArtifactGlLowTF, upRef, altAllele));

                VariableAllelicRatioGenotypeLikelihoods normalArtifactGlNF = normalReadPile.calculateLikelihoods(candidate.getNormalF(), normalReadPile.qualityScoreFilteredPileup);
                candidate.setNormalArtifactLodNF(normalReadPile.getAltVsRef(normalArtifactGlNF, upRef, altAllele));

                candidate.setInitialNormalAltQualitySum(normQs.getQualitySum(altAllele));
                candidate.setInitialNormalRefQualitySum(normQs.getQualitySum(upRef));

                candidate.setInitialNormalAltCounts(normQs.getCounts(altAllele));
                candidate.setInitialNormalRefCounts(normQs.getCounts(upRef));
                candidate.setInitialNormalReadDepth(normalReadPile.finalPileupReads.size());

                // TODO: don't hardcode.  Make this the max read length in the pile
                final long refStart = Math.max(1, rawContext.getLocation().getStart() - 150);
                String refGATKString = new String(ref.getBases());

                // handle the case where there are no bases before or after!
                char priorBasePositiveDirection = 'N';
                char priorBaseNegativeDirection = 'N';


                if (containsPosition(ref.getWindow(), candidate.getLocation().getStart() - 1)) {
                    priorBasePositiveDirection = (char) ref.getBases()[(int)candidate.getLocation().getStart() - MTAC.SEQ_ERROR_MODEL.getPriorBaseOffset() - (int)ref.getWindow().getStart()];
                }

                if (containsPosition(ref.getWindow(), candidate.getLocation().getStart() + 1)) {
                    priorBaseNegativeDirection = (char) ref.getBases()[(int)candidate.getLocation().getStart() + MTAC.SEQ_ERROR_MODEL.getPriorBaseOffset() - (int)ref.getWindow().getStart()];
                }

                candidate.setPriorBasePositiveDirection(priorBasePositiveDirection);
                candidate.setPriorBaseNegativeDirection(priorBaseNegativeDirection);

                // TODO: make this parameterizable
                final LocusReadPile t2 = filterReads(ref, tumorReadPile, true, !NO_BAQ);

                // if there are no reads remaining, abandon this theory
                if ( !MTAC.FORCE_OUTPUT && t2.finalPileupReads.size() == 0) { continue; }

                candidate.setInitialTumorAltCounts(t2.qualitySums.getCounts(altAllele));
                candidate.setInitialTumorRefCounts(t2.qualitySums.getCounts(upRef));
                candidate.setInitialTumorAltQualitySum(t2.qualitySums.getQualitySum(altAllele));
                candidate.setInitialTumorRefQualitySum(t2.qualitySums.getQualitySum(upRef));

                VariableAllelicRatioGenotypeLikelihoods t2Gl = t2.calculateLikelihoods(t2.finalPileup);
                candidate.setInitialTumorLod(t2.getAltVsRef(t2Gl, upRef, altAllele));
                candidate.setInitialTumorReadDepth(t2.finalPileupReads.size());

                candidate.setTumorF(t2.estimateAlleleFraction(upRef, altAllele));
                double tumorLod2 = t2.calculateAltVsRefLOD((byte)altAllele, candidate.getTumorF(), 0);
                candidate.setTumorLodFStar(tumorLod2);

                // calculate Tumor LOD with the local quality score
                candidate.setTumorLodLQS(t2.calculateAltVsRefLOD((byte)altAllele, candidate.getTumorF(), 0, lqs));

                // TODO: why extract the counts twice?  once above and once in this method...
                tumorFLB = 0;
                refCount = candidate.getInitialNormalRefCounts();
                altCount = candidate.getInitialTumorAltCounts();
                depth = refCount + altCount;
                if ( altCount > 0) {
                    // implemented as shown http://www.sigmazone.com/binomial_confidence_interval.htm
                    BetaDistribution dist = new BetaDistributionImpl(depth - altCount + 1, altCount);
                    tumorFLB = 1 - dist.inverseCumulativeProbability(1 - (1-0.95)/2);
                }
                candidate.setTumorFLowerBound(tumorFLB);

                //TODO: shouldn't this be f2 in the lod calculation instead of the strand specific f values?
                ReadBackedPileup forwardPileup = t2.finalPileup.getPositiveStrandPileup();
                double f2forward = LocusReadPile.estimateAlleleFraction(forwardPileup, upRef, altAllele);
                candidate.setTumorLodFStarForward(t2.calculateAltVsRefLOD(forwardPileup, (byte)altAllele, f2forward, 0.0, null));

                ReadBackedPileup reversePileup = t2.finalPileup.getNegativeStrandPileup();
                double f2reverse = LocusReadPile.estimateAlleleFraction(reversePileup, upRef, altAllele);
                candidate.setTumorLodFStarReverse(t2.calculateAltVsRefLOD(reversePileup, (byte)altAllele, f2reverse, 0.0, null));



                // calculate strand bias power
                candidate.setPowerToDetectPositiveStrandArtifact(
                        strandArtifactPowerCalculator.cachingPowerCalculation(reversePileup.depthOfCoverage(), candidate.getTumorF())
                );
                candidate.setPowerToDetectNegativeStrandArtifact(
                        strandArtifactPowerCalculator.cachingPowerCalculation(forwardPileup.depthOfCoverage(), candidate.getTumorF())
                );


                ArrayList<PileupElement> mutantPileupElements = new ArrayList<PileupElement>();
                ArrayList<PileupElement> referencePileupElements = new ArrayList<PileupElement>();


                for (PileupElement p : t2.finalPileup) {
                    final SAMRecord read = p.getRead();
                    final int offset = p.getOffset();

                    if (read.getReadString().charAt(offset) == altAllele) {
                        mutantPileupElements.add(p);
                    } else if (read.getReadString().charAt(offset) == upRef) {
                        referencePileupElements.add(p);
                    } else {
                        // just drop the read...
                    }
                }

                ReadBackedPileup mutantPileup =
                        new ReadBackedPileupImpl(rawContext.getLocation(), mutantPileupElements);

                ReadBackedPileup referencePileup =
                        new ReadBackedPileupImpl(rawContext.getLocation(), referencePileupElements);

                // TODO: shouldn't this be refAllele here?
                final LocusReadPile mutantPile = new LocusReadPile(mutantPileup, altAllele, 0, 0);
                final LocusReadPile refPile =  new LocusReadPile(referencePileup, altAllele, 0, 0);


                candidate.setPerfectStrandBias(calculatePerfectStrandBias(mutantPile));
                candidate.setStrandBias(calculateStrandBias(refPile, mutantPile));

                // start with just the tumor pile
                candidate.setTumorAltForwardOffsetsInRead(getForwardOffsetsInRead(mutantPileup));
                candidate.setTumorAltReverseOffsetsInRead(getReverseOffsetsInRead(mutantPileup));

                if (candidate.getTumorAltForwardOffsetsInRead().size() > 0) {
                    double[] offsets = convertIntegersToDoubles(candidate.getTumorAltForwardOffsetsInRead());
                    double median = getMedian(offsets);
                    candidate.setTumorForwardOffsetsInReadMedian(median);
                    candidate.setTumorForwardOffsetsInReadMad(calculateMAD(offsets, median));
                }


                if (candidate.getTumorAltReverseOffsetsInRead().size() > 0) {
                    double[] offsets = convertIntegersToDoubles(candidate.getTumorAltReverseOffsetsInRead());
                    double median = getMedian(offsets);
                    candidate.setTumorReverseOffsetsInReadMedian(median);
                    candidate.setTumorReverseOffsetsInReadMad(calculateMAD(offsets, median));
                }


                // test to see if the candidate should be rejected
                performRejection(candidate);

                String csOutput = callStatsGenerator.generateCallStats(candidate);
                if (MTAC.FORCE_ALLELES) {
                    out.println(csOutput);
                } else {
//                    System.out.println("putting in " + altAllele + " with lod " + candidate.getInitialTumorLod());
                    messageByTumorLod.put(candidate.getInitialTumorLod(), csOutput);
                }
            }

            // write out the call stats for the "best" candidate
            if (!messageByTumorLod.isEmpty()) {
                out.println(messageByTumorLod.lastEntry().getValue());
            }
            return -1;
        } catch (Throwable t) {
            System.err.println("Error processing " + rawContext.getContig() + ":" + rawContext.getPosition());
            t.printStackTrace(System.err);
            
            throw new RuntimeException(t);
        }
    }


    private String createSequenceContext(ReferenceContext ref, int size) {
        // create a context of 3 bases before, then 'x' then three bases after
        int offset = ref.getLocus().getStart() - ref.getWindow().getStart();
        StringBuilder sb = new StringBuilder(7);
        for(byte b : Arrays.copyOfRange(ref.getBases(), offset - size, offset)) {
            sb.append(Character.toUpperCase((char)b));
        }
        sb.append('x');
        for(byte b : Arrays.copyOfRange(ref.getBases(), offset + 1, offset + 1 + size)) {
            sb.append(Character.toUpperCase((char)b));
        }
        return sb.toString();
    }


    /**
     * Helper function that returns the phred-scaled base quality score we should use for calculating
     * likelihoods for a pileup element.  May return 0 to indicate that the observation is bad, and may
     * cap the quality score by the mapping quality of the read itself.
     *
     * @param p
     * @param ignoreBadBases
     * @param capBaseQualsAtMappingQual
     * @param minBaseQual
     * @return
     */
    private static byte qualToUse(PileupElement p, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, int minBaseQual) {
        if ( ignoreBadBases && !BaseUtils.isRegularBase(p.getBase()) )
            return 0;

        byte qual = p.getQual();

        if ( qual > SAMUtils.MAX_PHRED_SCORE )
            throw new UserException.MalformedBAM(p.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, p.getRead().getReadName()));
        if ( capBaseQualsAtMappingQual )
            qual = (byte)Math.min((int)p.getQual(), p.getMappingQual());
        if ( (int)qual < minBaseQual )
            qual = (byte)0;

        return qual;
    }

    private boolean containsPosition(GenomeLoc window, int position) {
        return (window.getStart() <= position && position <= window.getStop());
    }

    public static class FisherData {
        public int a;
        public int b;
        public int c;
        public int d;
        public double p;
        public double maxPosP;
        public double maxNegP;

        public FisherData(int a, int b, int c, int d, double p, double maxPosP, double maxNegP) {
            this.a = a;
            this.b = b;
            this.c = c;
            this.d = d;
            this.p = p;
            this.maxPosP = maxPosP;
            this.maxNegP = maxNegP;
        }

        public String dataToString() {
            StringBuilder sb = new StringBuilder();
            sb.append("(");
            sb.append(a).append(",");
            sb.append(b).append(",");
            sb.append(c).append(",");
            sb.append(d).append(")");
            return sb.toString();
        }

        public int getA() {
            return a;
        }

        public int getB() {
            return b;
        }

        public int getC() {
            return c;
        }

        public int getD() {
            return d;
        }

        public double getP() {
            return p;
        }

        public double getMaxPosP() {
            return maxPosP;
        }

        public void setMaxPosP(double maxPosP) {
            this.maxPosP = maxPosP;
        }

        public double getMaxNegP() {
            return maxNegP;
        }

        public void setMaxNegP(double maxNegP) {
            this.maxNegP = maxNegP;
        }
    }

    private FisherData calculateClippingBias(LocusReadPile refPile, LocusReadPile mutantPile) {
        // Construct a 2x2 contingency table of
        //            not-clipped  clipped
        //      REF    a            b
        //      MUT    c            d
        //
        // and return an array of {a,b,c,d,two-tailed-p}

        int a = 0, b = 0, c = 0, d = 0;
        for (SAMRecord rec : refPile.finalPileupReads) {
            if (isReadHeavilySoftClipped(rec)) { b++;} else { a++; }
        }
        for (SAMRecord rec : mutantPile.finalPileupReads) {
            if (isReadHeavilySoftClipped(rec)) { d++;} else { c++; }
        }

        final double p = fisher.getTwoTailedP(a,b,c,d);
        int n = c+d;
        final double pmaxpos = fisher.getTwoTailedP(a,b,n,0);
        final double pmaxneg = fisher.getTwoTailedP(a,b,0,n);

        return new FisherData(a,b,c,d,p, pmaxpos, pmaxneg);
    }

    // TODO: we can do this more cheaply with a GATKSAMRecord...
    private boolean isReadHeavilySoftClipped(SAMRecord rec) {
        int total = 0;
        int clipped = 0;
        for(CigarElement ce : rec.getCigar().getCigarElements()) {
            total += ce.getLength();
            if (ce.getOperator() == CigarOperator.SOFT_CLIP) {
                clipped += ce.getLength();
            }
        }

        return ((float) clipped / (float)total >= MTAC.HEAVILY_CLIPPED_READ_FRACTION);
    }

    private FisherData calculateStrandBias(LocusReadPile refPile, LocusReadPile mutantPile) {
        // Construct a 2x2 contingency table of
        //            pos     neg
        //      REF    a       b
        //      MUT    c       d
        //
        // and return an array of {a,b,c,d,two-tailed-p}

        //
        // POWER: to meet P <= 0.05 with just 1 mut allele we need 19 ref reads (f=0.05)
        // to meet P <= 0.05 with 2 mut alleles, we need 5 ref reads
        // to meet P <= 0.05 with 3 mut alleles, we need 3 ref reads 

        
        int a = 0, b = 0, c = 0, d = 0;
        for (SAMRecord rec : refPile.finalPileupReads) {
            if (rec.getReadNegativeStrandFlag()) { b++;} else { a++; }
        }
        for (SAMRecord rec : mutantPile.finalPileupReads) {
            if (rec.getReadNegativeStrandFlag()) { d++;} else { c++; }
        }

        final double p = fisher.getTwoTailedP(a,b,c,d);
        int n = c+d;
        final double pmaxpos = fisher.getTwoTailedP(a,b,n,0);
        final double pmaxneg = fisher.getTwoTailedP(a,b,0,n);

        return new FisherData(a,b,c,d,p, pmaxpos, pmaxneg);
    }


    private FisherData calculatePerfectStrandBias(LocusReadPile mutantPile) {
        try {
            int c = 0, d = 0;
            for (SAMRecord rec : mutantPile.finalPileupReads) {
                if (rec.getReadNegativeStrandFlag()) { d++;} else { c++; }
            }

            int trials = c + d;
            BinomialDistribution bd = new BinomialDistributionImpl(trials, 0.5d);

            double pmax = bd.cumulativeProbability(0);
            double p = bd.cumulativeProbability(Math.min(c,d));

            return new FisherData(0,0,c,d,p,pmax, pmax);
        } catch (MathException me) {
            throw new RuntimeException("Error in ApacheMath"+me.getMessage(), me);
        }
    }

    private List<Integer> getForwardOffsetsInRead(ReadBackedPileup p) {
        return getOffsetsInRead(p, true);
    }

    private List<Integer> getReverseOffsetsInRead(ReadBackedPileup p) {
        return getOffsetsInRead(p, false);
    }

    private List<Integer> getOffsetsInRead(ReadBackedPileup p, boolean useForwardOffsets) {
        List<Integer> positions = new ArrayList<Integer>();
        for(PileupElement pe : p) {
            // TODO: maybe we should be doing start-site distribution, or a clipping aware offset?
                GATKSAMRecord r = pe.getRead();

                positions.add(
                        Math.abs((int)(p.getLocation().getStart() -
                                (useForwardOffsets?r.getOriginalAlignmentStart():r.getOriginalAlignmentEnd())))
                );
        }

        return positions;
    }

    private static final double PERFECT_STRAND_BIAS_THRESHOLD = .05d;
    private static final double STRAND_BIAS_THRESHOLD = .05d;

    private static final int QUALITY_RST_MEDIAN_SHIFT_THRESHOLD = 4;
    private static final double QUALITY_RST_PVALUE_THRESHOLD = .005;

    private void performRejection(CandidateMutation candidate) {
        if (candidate.getTumorLodFStar() < MTAC.TUMOR_LOD_THRESHOLD) {
            candidate.addRejectionReason("fstar_tumor_lod");
        }

        if (MTAC.ARTIFACT_DETECTION_MODE) {
            return;
        }

        // if the best theory for the normal is A het (not necessarily this het)
        // just move on.  We're not attempting to call LOH with this tool
        if (candidate.getInitialNormalBestGenotype() != DiploidGenotype.AA &&
            candidate.getInitialNormalBestGenotype() != DiploidGenotype.CC &&
            candidate.getInitialNormalBestGenotype() != DiploidGenotype.GG &&
            candidate.getInitialNormalBestGenotype() != DiploidGenotype.TT) {
            candidate.addRejectionReason("het_normal");
        }

        if (candidate.getTumorInsertionCount() >= MTAC.GAP_EVENTS_THRESHOLD ||
            candidate.getTumorDeletionCount()  >= MTAC.GAP_EVENTS_THRESHOLD) {
            candidate.addRejectionReason("nearby_gap_events");
        }

        if (MTAC.FRACTION_CONTAMINATION+MTAC.MINIMUM_MUTATION_CELL_FRACTION > 0 && candidate.getTumorLodFStar() <= MTAC.TUMOR_LOD_THRESHOLD + Math.max(0, candidate.getContaminantLod())) {
            candidate.addRejectionReason("possible_contamination");
        }

        //TODO: clean up this rejection reason. no space and it's really "germline risk" or something
        if (candidate.isGermlineAtRisk() && candidate.getInitialNormalLod() < MTAC.NORMAL_DBSNP_LOD_THRESHOLD) {
            candidate.addRejectionReason("DBSNP Site");
        }

        if (candidate.getInitialNormalLod() < MTAC.NORMAL_LOD_THRESHOLD) {
            candidate.addRejectionReason("normal_lod");
        }

//        if ( (candidate.getInitialNormalAltCounts() >= MAX_ALT_ALLELES_IN_NORMAL_COUNT && candidate.getInitialNormalAltQualitySum() > MAX_ALT_ALLELES_IN_NORMAL_QSCORE_SUM && candidate.getNormalF() > MAX_ALT_ALLELE_IN_NORMAL_FRACTION)) {
//            candidate.addRejectionReason("alt allele in normal");
//        }

        if (candidate.getNormalArtifactLodNF() > MTAC.NORMAL_ARTIFACT_LOD_THRESHOLD) {
            candidate.addRejectionReason("normal_artifact_lod");
        }

        if ( (candidate.getTumorForwardOffsetsInReadMedian() != null && candidate.getTumorForwardOffsetsInReadMedian() <= MTAC.PIR_MEDIAN_THRESHOLD && candidate.getTumorForwardOffsetsInReadMad() != null && candidate.getTumorForwardOffsetsInReadMad() <= MTAC.PIR_MAD_THRESHOLD) ||
              candidate.getTumorReverseOffsetsInReadMedian() != null && candidate.getTumorReverseOffsetsInReadMedian() <= MTAC.PIR_MEDIAN_THRESHOLD && candidate.getTumorReverseOffsetsInReadMad() != null && candidate.getTumorReverseOffsetsInReadMad() <= MTAC.PIR_MAD_THRESHOLD ) {
            candidate.addRejectionReason("clustered_read_position");

        }


        candidate.setPositiveDirectionPowered(
                candidate.getPerfectStrandBias().getMaxPosP() < PERFECT_STRAND_BIAS_THRESHOLD &&
                candidate.getStrandBias().getMaxPosP() < STRAND_BIAS_THRESHOLD);

        candidate.setNegativeDirectionPowered(
                candidate.getPerfectStrandBias().getMaxNegP() < PERFECT_STRAND_BIAS_THRESHOLD &&
                candidate.getStrandBias().getMaxNegP() < STRAND_BIAS_THRESHOLD);


        if (candidate.getPerfectStrandBias().getP() == Double.NaN ||
            candidate.getStrandBias().getP() == Double.NaN) {
            candidate.addRejectionReason("ERROR: Unable to calculate Strand Bias Score");
        }

        // TODO: parameterize this 2.0 value
        // TODO: sync naming (is it positive or forward)?
        // test both the positive and negative direction.  If you're "at risk" for an artifact in a direction you need
        // to observerve a LOD score > than the threshold in the OPPOSITE direction
//        if (candidate.isPositiveDirectionAtRisk() && candidate.isPositiveDirectionPowered() && candidate.getTumorLodFStarReverse() < 2.0) {
//            candidate.addRejectionReason("positive_strand_artifact");
//        }
//
//        if (candidate.isNegativeDirectionAtRisk() && candidate.isNegativeDirectionPowered() && candidate.getTumorLodFStarForward() < 2.0) {
//            candidate.addRejectionReason("negative_strand_artifact");
//        }

        if (candidate.getPowerToDetectNegativeStrandArtifact() >= MTAC.STRAND_ARTIFACT_POWER_THRESHOLD && candidate.getTumorLodFStarForward() < MTAC.STRAND_ARTIFACT_LOD_THRESHOLD) {
            candidate.addRejectionReason("strand_artifact");
        }

        if (candidate.getPowerToDetectPositiveStrandArtifact() >= MTAC.STRAND_ARTIFACT_POWER_THRESHOLD && candidate.getTumorLodFStarReverse() < MTAC.STRAND_ARTIFACT_LOD_THRESHOLD) {
            candidate.addRejectionReason("strand_artifact");
        }

        if (candidate.getTotalPairs() > 0 && ((float)candidate.getMapQ0Reads() / (float)candidate.getTotalPairs()) >= MTAC.FRACTION_MAPQ0_THRESHOLD) {
            candidate.addRejectionReason("poor_mapping_region_mapq0");
        }

        if (candidate.isSeenInPanelOfNormals()) {
            if (candidate.isCosmicSite()) {
                // if we saw it in the panel of normals, retain the call it was a COSMIC, but non-dbsnp site,
            } else {
                // otherwise, reject it
                candidate.addRejectionReason("seen_in_panel_of_normals");
            }
        }
    }

    private double calculateMAD(double[] dd, double median) {
        double[] dev = new double[dd.length];
        for(int i=0; i<dd.length; i++) {
            dev[i] = Math.abs(dd[i] - median);
        }
        return getMedian(dev);

    }

    private double getMedian(double[] data) {
        Arrays.sort(data);
        Double result;

        if (data.length % 2 == 1) {
            // If the number of entries in the list is not even.

            // Get the middle value.
            // You must floor the result of the division to drop the
            // remainder.
            result = data[(int) Math.floor(data.length/2) ];

        } else {
            // If the number of entries in the list are even.

            // Get the middle two values and average them.
            Double lowerMiddle = data[data.length/2 ];
            Double upperMiddle = data[data.length/2 - 1 ];
            result = (lowerMiddle + upperMiddle) / 2;
        }

        return result;
    }

    public static double[] convertIntegersToDoubles(List<Integer> integers)
    {
        double[] ret = new double[integers.size()];
        for (int i=0; i < ret.length; i++)
        {
            ret[i] = integers.get(i);
        }
        return ret;
    }

    public Integer treeReduce(Integer lhs, Integer rhs) {
        return 0;
    }

    // Given result of map function
    @Override
	public Integer reduceInit() {
        return 0;
    }
    @Override
	public Integer reduce(final Integer value, final Integer sum) {
        return 0;
    }



    int MAX_READ_MISMATCH_QUALITY_SCORE_SUM = 100;
    private static Character MAPPED_BY_MATE = 'M';
    private BAQ baqHMM = new BAQ();
    // TODO: or should we be using this? baqHMM = new BAQ(1e-3, 0.1, bw, (byte)0, true); from the BAQ unit test?
    IndexedFastaSequenceFile refReader;

    private LocusReadPile filterReads(final ReferenceContext ref, final LocusReadPile pile, boolean filterMateRescueReads, boolean applyBAQ) {
        ArrayList<PileupElement> newPileupElements = new ArrayList<PileupElement>();

        for ( PileupElement p : pile.finalPileup ) {
            final GATKSAMRecord read = p.getRead();

            int mismatchQualitySum =
                    AlignmentUtils.mismatchesInRefWindow(p, ref, false, true);

            // do we have to many mismatches overall?
            if (mismatchQualitySum > this.MAX_READ_MISMATCH_QUALITY_SCORE_SUM) {
                //if (read.getReadString().charAt(offset) == altAllele) { out.println("MAXERR " + read.getReadName() + " that supported altAllele with " + mismatchScore + " mmqsum"); }
                continue;
            }

            // is this a heavily clipped read?
            if (isReadHeavilySoftClipped(read)) {
                continue;
            }

            // was this read ONLY placed because it's mate was uniquely placed? (supplied by BWA)
            if (filterMateRescueReads && MAPPED_BY_MATE.equals(read.getAttribute("XT"))) {
                continue;
            }

            // apply BAQ
            if (applyBAQ) {
                baqHMM.baqRead(read, refReader, BAQ.CalculationMode.RECALCULATE, BAQ.QualityMode.OVERWRITE_QUALS);
            }

            // if we're here... we passed all the read filters!
            newPileupElements.add(new PileupElement(read, p.getOffset(), p.isDeletion(), p.isBeforeDeletionStart(), p.isAfterDeletionEnd(), p.isBeforeInsertion(), p.isAfterInsertion(),p.isNextToSoftClip()));


        }
        ReadBackedPileup newPileup =
                new ReadBackedPileupImpl(pile.pileup.getLocation(), newPileupElements);


        final LocusReadPile newPile = new LocusReadPile(newPileup, pile.refBase, 0, 0);

        return newPile;
    }

    /**
     * Calculate the fraction of mismatches to the altAllele against all other
     * alleles from the offset position in the read to the end of the read
     *
     * @param altAllele
     * @return
     */
    private float getAltAlleleMismatchRate(final SAMRecord read, final List<Mismatch> mismatches, final int offset, final char altAllele) {
        int mmToAlt = 0;
        int mmToOther = 0;
        for (final Mismatch mm : mismatches) {
            if ((read.getReadNegativeStrandFlag() && mm.offset < offset) ||
                (!read.getReadNegativeStrandFlag() && mm.offset > offset) ) {
                if (mm.mismatchBase == altAllele) {
                    mmToAlt++;
                } else {
                    mmToOther++;
                }
            }
        }

        // TODO: use fischer exact test?
        
        return (float) (mmToAlt + 1) / (float) (mmToAlt + mmToOther + 4);
    }





    // ----------------------------------- PRIVATE IN IntervalCleanerWalker
    public static final int MAX_QUAL = 99;

    private static class Mismatch {
        SAMRecord read;
        int offset;
        long position;
        char mismatchBase;
        int qualityScore;

        private Mismatch(final SAMRecord read, final int offset, final long position, final char mismatchBase, final int qualityScore) {
            this.read = read;
            this.offset = offset;
            this.position = position;
            this.mismatchBase = mismatchBase;
            this.qualityScore = qualityScore;
        }
    }

    private int mismatchQualitySum(final SAMRecord read, final String refSeq, final int refIndex, final int minMismatchQualityScore) {
        final List<Mismatch> mismatches = getMismatches(read, refSeq, refIndex);
        return mismatchQualitySum(mismatches, minMismatchQualityScore);
    }

    private int mismatchQualitySum(final List<Mismatch> mismatches, final int minMismatchQualityScore) {
        int sum = 0;
        for(final Mismatch mm : mismatches) {
            if (mm.qualityScore >= minMismatchQualityScore) {
                sum += mm.qualityScore;
            }
        }
        return sum;
    }




    /**
     * Returns a list of position
     * @param refSeq
     * @param refIndex
     * @return
     */
    private static List<Mismatch> getMismatches(final SAMRecord read, final String refSeq, int refIndex) {
        final List<Mismatch> mismatches = new ArrayList<Mismatch>();

        final String readSeq = read.getReadString();
        final String quals = read.getBaseQualityString();
        int readIndex = 0;
        int sum = 0;
        final Cigar c = read.getCigar();
        for (int i = 0 ; i < c.numCigarElements() ; i++) {
            final CigarElement ce = c.getCigarElement(i);
            switch ( ce.getOperator() ) {
                case M:
                    for (int j = 0 ; j < ce.getLength() ; j++, refIndex++, readIndex++ ) {
                        // TODO: what is this case????
                        if ( refIndex >= refSeq.length() ) {
                            sum += MAX_QUAL;
                        } else {
                            final char readBase = Character.toUpperCase(readSeq.charAt(readIndex));
                            final char refBase = Character.toUpperCase(refSeq.charAt(refIndex));

                            if ( readBase != refBase ) {
                                final int qual = quals.charAt(readIndex) - 33;
                                mismatches.add(new Mismatch(read, readIndex, refIndex, readBase, qual));
                            }
                        }
                    }
                    break;
                case I:
                    readIndex += ce.getLength();
                    break;
                case D:
                    refIndex += ce.getLength();
                    break;
            }

        }
        return mismatches;
    }

    static LinkedHashMap sortByDescendingValue(Map map) {
     List list = new LinkedList(map.entrySet());
     Collections.sort(list, new Comparator() {
          public int compare(Object o1, Object o2) {
               return -1 * ((Comparable) ((Map.Entry) (o1)).getValue())
              .compareTo(((Map.Entry) (o2)).getValue());
          }
     });
    // logger.info(list);
    LinkedHashMap result = new LinkedHashMap();
    for (Iterator it = list.iterator(); it.hasNext();) {
        Map.Entry entry = (Map.Entry)it.next();
        result.put(entry.getKey(), entry.getValue());
     }
    return result;
}

    private boolean isTumorRead(SAMRecord read) {
        SAMReaderID id = getToolkit().getReaderIDForRead(read);
        return tumorSAMReaderIDs.contains(id);
    }

}
