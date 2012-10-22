package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import java.io.*;
import java.util.*;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.*;

import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.variantcontext.*;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;
import org.apache.commons.lang.math.NumberUtils;
import org.apache.commons.math.distribution.BetaDistribution;
import org.apache.commons.math.distribution.BetaDistributionImpl;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

@PartitionBy(PartitionType.LOCUS)
@BAQMode()
@Reference(window=@Window(start=-1* MuTect.REFERENCE_HALF_WINDOW_LENGTH,stop= MuTect.REFERENCE_HALF_WINDOW_LENGTH))
@By(DataSource.REFERENCE)
public class MuTect extends LocusWalker<Integer, Integer> implements TreeReducible<Integer> {
    public static final int REFERENCE_HALF_WINDOW_LENGTH = 150;
    public static final String BAM_TAG_TUMOR = "tumor";
    public static final String BAM_TAG_NORMAL = "normal";
    public static final String BAM_TAG_CONTROL = "control";

    // DO NOT CHANGE THIS LINE!  It's the SVN revision number of the caller, which updates automatically!
    private static final String VERSION = "$Rev$";

    @ArgumentCollection private MuTectArgumentCollection MTAC = new MuTectArgumentCollection();

    /***************************************/
    // Call-stats output
    /***************************************/
    @Output(doc="Call-stats output")
    PrintStream out;

    @Output(doc="VCF output of mutation candidates",shortName="vcf", fullName="vcf", required=false)
    protected VariantContextWriter vcf = null;

    /***************************************/
    // Reference Metadata inputs
    /***************************************/
    @Input(fullName="dbsnp", shortName = "dbsnp", doc="VCF file of DBSNP information", required=false)
    public List<RodBinding<VariantContext>> dbsnpRod = Collections.emptyList();;

    @Input(fullName="cosmic", shortName = "cosmic", doc="VCF file of COSMIC sites", required=false)
    public List<RodBinding<VariantContext>> cosmicRod = Collections.emptyList();;

    @Hidden
    @Input(fullName="normal_panel", shortName = "normal_panel", doc="VCF file of sites observed in normal", required=false)
    public List<RodBinding<VariantContext>> normalPanelRod = Collections.emptyList();;

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


    private CoverageWiggleFileWriter stdCovWriter;
    private CoverageWiggleFileWriter q20CovWriter;
    private CoverageWiggleFileWriter powerWriter;
    private CoverageWiggleFileWriter tumorDepthWriter;
    private CoverageWiggleFileWriter normalDepthWriter;


    private Set<SAMReaderID> tumorSAMReaderIDs = new HashSet<SAMReaderID>();
    private Set<SAMReaderID> normalSAMReaderIDs = new HashSet<SAMReaderID>();
    private Set<SAMReaderID> controlSAMReaderIDs = new HashSet<SAMReaderID>();

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
                throw new RuntimeException("BAMs must be tagged as either 'tumor','normal' or 'control'");
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
                    normalSAMReaderIDs.add(id);


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
                } else if (BAM_TAG_CONTROL.equalsIgnoreCase(tag)) {                    
                    controlSAMReaderIDs.add(id);
                } else {
                    throw new RuntimeException("Unknown BAM tag '" + tag + "' must be either 'tumor','normal' or 'control'");
                }                
            }
        }

        if (!hasTumorBam) {
            throw new RuntimeException("At least one BAM tagged as 'tumor' required");
        }

        // FIXME: update for Control BAM concept
        if (!hasNormalBam) {
            MTAC.NORMAL_LOD_THRESHOLD = -1 * Float.MAX_VALUE;
            MTAC.NORMAL_DBSNP_LOD_THRESHOLD = -1 * Float.MAX_VALUE;
            MTAC.NORMAL_ARTIFACT_LOD_THRESHOLD = Float.MAX_VALUE;
            MTAC.NORMAL_SAMPLE_NAME = "none";
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

        // initialize the call-stats file
        out.println("## muTector v1.0." + VERSION.split(" ")[1]);
        out.println(callStatsGenerator.generateHeader());

        // initialize the VCF output
        if (vcf != null) {
            Set<String> samples = new HashSet<String>();
            samples.add(MTAC.TUMOR_SAMPLE_NAME);
            samples.add(MTAC.NORMAL_SAMPLE_NAME);
            Set<VCFHeaderLine> headerInfo = getVCFHeaderInfo(MTAC);
            vcf.writeHeader(new VCFHeader(headerInfo, samples));
        }

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

        TreeMap<Double, CandidateMutation> messageByTumorLod = new TreeMap<Double, CandidateMutation>();

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

            ArrayList<PileupElement> tumorPileupElements = new ArrayList<PileupElement>();
            ArrayList<PileupElement> normalPileupElements = new ArrayList<PileupElement>();
            ArrayList<PileupElement> controlPileupElements = new ArrayList<PileupElement>();

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
                ReadSource source = getReadSource(p.getRead());
                if (source == ReadSource.Tumor) {
                    tumorPileupElements.add(p);
                } else if (source == ReadSource.Normal) {
                    normalPileupElements.add(p);
                } else if (source == ReadSource.Control) {
                    controlPileupElements.add(p);
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
                candidate.setPanelOfNormalsVC(panelOfNormalsVC.isEmpty()?null:panelOfNormalsVC.iterator().next()); // if there are multiple, we're just grabbing the first
                candidate.setCosmicSite(!cosmicVC.isEmpty());
                candidate.setDbsnpSite(!dbsnpVC.isEmpty());
                candidate.setDbsnpVC(dbsnpVC.isEmpty()?null:dbsnpVC.iterator().next());
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

                candidate.setNormalFQuals(LocusReadPile.estimateAlleleFractionUsingQuals(normalReadPile.qualityScoreFilteredPileup, (byte)upRef, (byte)altAllele));
                VariableAllelicRatioGenotypeLikelihoods normalArtifactGlNFQ = normalReadPile.calculateLikelihoods(candidate.getNormalFQuals(), normalReadPile.qualityScoreFilteredPileup);
                candidate.setNormalArtifactLodNFQ(normalReadPile.getAltVsRef(normalArtifactGlNFQ, upRef, altAllele));


                candidate.setInitialNormalAltQualitySum(normQs.getQualitySum(altAllele));
                candidate.setInitialNormalRefQualitySum(normQs.getQualitySum(upRef));

                candidate.setInitialNormalAltCounts(normQs.getCounts(altAllele));
                candidate.setInitialNormalRefCounts(normQs.getCounts(upRef));
                candidate.setInitialNormalReadDepth(normalReadPile.finalPileupReads.size());

                // TODO: make this parameterizable
                final LocusReadPile t2 = filterReads(ref, tumorReadPile.finalPileup, true);

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
                refCount = candidate.getInitialTumorRefCounts();
                altCount = candidate.getInitialTumorAltCounts();
                depth = refCount + altCount;
                if ( altCount > 0) {
                    // implemented as shown http://www.sigmazone.com/binomial_confidence_interval.htm
                    BetaDistribution dist = new BetaDistributionImpl(depth - altCount + 1, altCount);
                    tumorFLB = 1 - dist.inverseCumulativeProbability(1 - (1-0.95)/2);
                }
                candidate.setTumorFLowerBound(tumorFLB);

                //TODO: shouldn't this be f2 in the lod calculation instead of the strand specific f values?
                // TODO: clean up use of forward/reverse vs positive/negative (prefer the latter since GATK uses it)
                ReadBackedPileup forwardPileup = filterReads(ref, tumorReadPile.finalPileupPositiveStrand, true).finalPileupPositiveStrand;
                double f2forward = LocusReadPile.estimateAlleleFraction(forwardPileup, upRef, altAllele);
                candidate.setTumorLodFStarForward(t2.calculateAltVsRefLOD(forwardPileup, (byte)altAllele, f2forward, 0.0, null));

                ReadBackedPileup reversePileup = filterReads(ref,tumorReadPile.finalPileupNegativeStrand, true).finalPileupNegativeStrand;
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

                // Set the maximum observed mapping quality score for the reference and alternate alleles
                byte[] rmq = referencePileup.getMappingQuals();
                candidate.setTumorRefMaxMapQ((rmq.length==0)?0:NumberUtils.max(rmq));

                byte[] amq = mutantPileup.getMappingQuals();
                candidate.setTumorAltMaxMapQ((amq.length==0)?0:NumberUtils.max(amq));

                candidate.setStrandContingencyTable(getStrandContingencyTable(refPile, mutantPile));

                // start with just the tumor pile
                candidate.setTumorAltForwardOffsetsInRead(getForwardOffsetsInRead(mutantPileup));
                candidate.setTumorAltReverseOffsetsInRead(getReverseOffsetsInRead(mutantPileup));

                if (candidate.getTumorAltForwardOffsetsInRead().size() > 0) {
                    double[] offsets = MuTectStats.convertIntegersToDoubles(candidate.getTumorAltForwardOffsetsInRead());
                    double median = MuTectStats.getMedian(offsets);
                    candidate.setTumorForwardOffsetsInReadMedian(median);
                    candidate.setTumorForwardOffsetsInReadMad(MuTectStats.calculateMAD(offsets, median));
                }


                if (candidate.getTumorAltReverseOffsetsInRead().size() > 0) {
                    double[] offsets = MuTectStats.convertIntegersToDoubles(candidate.getTumorAltReverseOffsetsInRead());
                    double median = MuTectStats.getMedian(offsets);
                    candidate.setTumorReverseOffsetsInReadMedian(median);
                    candidate.setTumorReverseOffsetsInReadMad(MuTectStats.calculateMAD(offsets, median));
                }


                // test to see if the candidate should be rejected
                performRejection(candidate);

                // FIXME: this is inefficient.  Put everything into the data structure (with Candidate as the value) and then print it outside the main loop
                if (MTAC.FORCE_ALLELES) {
                    out.println(callStatsGenerator.generateCallStats(candidate));
                    // TODO: should force alleles be enabled in VCF
                } else {
                    messageByTumorLod.put(candidate.getInitialTumorLod(), candidate);
                }
            }

            // if more than one site passes the tumor lod threshold for KEEP the fail the tri_allelic Site filter
            int passingCandidates = 0;
            for(CandidateMutation c : messageByTumorLod.values()) {
                if (c.getTumorLodFStar() >= MTAC.TUMOR_LOD_THRESHOLD){
                    passingCandidates++;
                }
            }

            if (passingCandidates > 1) {
                for(CandidateMutation c : messageByTumorLod.values()) {
                    c.addRejectionReason("triallelic_site");
                }
            }

            // write out the call stats for the "best" candidate
            if (!messageByTumorLod.isEmpty()) {
                CandidateMutation m = messageByTumorLod.lastEntry().getValue();

                // only output passing calls OR rejected sites if ONLY_PASSING_CALLS is not specified
                if (!m.isRejected() || (m.isRejected() && !MTAC.ONLY_PASSING_CALLS)) {

                    out.println(callStatsGenerator.generateCallStats(m));
                    if (vcf != null) {
                        VariantContext vc = generateVC(m);
                        vcf.add(vc);
                    }
                }
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

        for(byte b : Arrays.copyOfRange(ref.getBases(), Math.max(0,offset - size), offset)) {
            sb.append(Character.toUpperCase((char)b));
        }
        sb.append('x');
        for(byte b : Arrays.copyOfRange(ref.getBases(), offset + 1, Math.min(ref.getBases().length,offset + 1 + size))) {
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

    private int[] getStrandContingencyTable(LocusReadPile refPile, LocusReadPile mutantPile) {
        // Construct a 2x2 contingency table of
        //            pos     neg
        //      REF    a       b
        //      MUT    c       d
        //
        // and return an array of {a,b,c,d}

        int a = 0, b = 0, c = 0, d = 0;
        for (SAMRecord rec : refPile.finalPileupReads) {
            if (rec.getReadNegativeStrandFlag()) { b++;} else { a++; }
        }
        for (SAMRecord rec : mutantPile.finalPileupReads) {
            if (rec.getReadNegativeStrandFlag()) { d++;} else { c++; }
        }

        int[] results = new int[]{a,b,c,d};
        return results;
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
                        Math.abs((int)(p.getLocation().getStart() - (useForwardOffsets?pe.getRead().getAlignmentStart():pe.getRead().getAlignmentEnd())))
// TODO: the following code was for handling reduced reads, but if you use it on non-reduced reads it returns the wrong thing AND there is no way to tell if this is a reduced read!
//                        Math.abs((int)(p.getLocation().getStart() -
//                                (useForwardOffsets?r.getOriginalAlignmentStart():r.getOriginalAlignmentEnd())))
                );
        }

        return positions;
    }

    private void performRejection(CandidateMutation candidate) {
        if (candidate.getTumorLodFStar() < MTAC.TUMOR_LOD_THRESHOLD) {
            candidate.addRejectionReason("fstar_tumor_lod");
        }

        if (MTAC.ARTIFACT_DETECTION_MODE) {
            return;
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

        if ( (candidate.getInitialNormalAltCounts() >= MTAC.MAX_ALT_ALLELES_IN_NORMAL_COUNT || candidate.getNormalF() >= MTAC.MAX_ALT_ALLELE_IN_NORMAL_FRACTION ) && candidate.getInitialNormalAltQualitySum() > MTAC.MAX_ALT_ALLELES_IN_NORMAL_QSCORE_SUM) {
            candidate.addRejectionReason("alt_allele_in_normal");
        }

        if ( (candidate.getTumorForwardOffsetsInReadMedian() != null && candidate.getTumorForwardOffsetsInReadMedian() <= MTAC.PIR_MEDIAN_THRESHOLD && candidate.getTumorForwardOffsetsInReadMad() != null && candidate.getTumorForwardOffsetsInReadMad() <= MTAC.PIR_MAD_THRESHOLD) ||
              candidate.getTumorReverseOffsetsInReadMedian() != null && candidate.getTumorReverseOffsetsInReadMedian() <= MTAC.PIR_MEDIAN_THRESHOLD && candidate.getTumorReverseOffsetsInReadMad() != null && candidate.getTumorReverseOffsetsInReadMad() <= MTAC.PIR_MAD_THRESHOLD ) {
            candidate.addRejectionReason("clustered_read_position");

        }

        // TODO: sync naming (is it positive or forward)?
        if (
                (candidate.getPowerToDetectNegativeStrandArtifact() >= MTAC.STRAND_ARTIFACT_POWER_THRESHOLD && candidate.getTumorLodFStarForward() < MTAC.STRAND_ARTIFACT_LOD_THRESHOLD) ||
                (candidate.getPowerToDetectPositiveStrandArtifact() >= MTAC.STRAND_ARTIFACT_POWER_THRESHOLD && candidate.getTumorLodFStarReverse() < MTAC.STRAND_ARTIFACT_LOD_THRESHOLD)
                ) {
            candidate.addRejectionReason("strand_artifact");
        }

        if (candidate.getTotalPairs() > 0 && ((float)candidate.getMapQ0Reads() / (float)candidate.getTotalPairs()) >= MTAC.FRACTION_MAPQ0_THRESHOLD) {
            candidate.addRejectionReason("poor_mapping_region_mapq0");
        }
        
        if (candidate.getTumorAltMaxMapQ() < MTAC.REQUIRED_MAXIMUM_ALT_ALLELE_MAPPING_QUALITY_SCORE) {
            candidate.addRejectionReason("poor_mapping_region_alternate_allele_mapq");
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
    // TODO: or should we be using this? baqHMM = new BAQ(1e-3, 0.1, bw, (byte)0, true); from the BAQ unit test?
    IndexedFastaSequenceFile refReader;

    private LocusReadPile filterReads(final ReferenceContext ref, final ReadBackedPileup pile, boolean filterMateRescueReads) {
        ArrayList<PileupElement> newPileupElements = new ArrayList<PileupElement>();

        for ( PileupElement p : pile ) {
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

            // if we're here... we passed all the read filters!
            newPileupElements.add(new PileupElement(read, p.getOffset(), p.isDeletion(), p.isBeforeDeletionStart(), p.isAfterDeletionEnd(), p.isBeforeInsertion(), p.isAfterInsertion(),p.isNextToSoftClip()));


        }
        ReadBackedPileup newPileup =
                new ReadBackedPileupImpl(ref.getLocus(), newPileupElements);


        final LocusReadPile newPile = new LocusReadPile(newPileup, (char)ref.getBase(), 0, 0);

        return newPile;
    }

    public enum ReadSource { Tumor, Normal, Control }
    
    private ReadSource getReadSource(SAMRecord read) {
        // check if it's a tumor
        SAMReaderID id = getToolkit().getReaderIDForRead(read);
        if (tumorSAMReaderIDs.contains(id)) { return ReadSource.Tumor; }
        if (normalSAMReaderIDs.contains(id)) { return ReadSource.Normal; }
        if (controlSAMReaderIDs.contains(id)) { return ReadSource.Control; }
        
        // unexpected condition
        throw new RuntimeException("Unable to determine read source (tumor,normal,control) for read " + read.getReadName());               
    }

    private static Set<VCFHeaderLine> getVCFHeaderInfo(final MuTectArgumentCollection MTAC) {
        Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();


        headerInfo.add(new VCFFilterHeaderLine("REJECT", "Rejected as a confident somatic mutation"));
        headerInfo.add(new VCFFilterHeaderLine("PASS", "Accept as a confident somatic mutation"));

        // TODO: what fields do we need here
        VCFStandardHeaderLines.addStandardInfoLines(headerInfo, true,
                VCFConstants.MAPPING_QUALITY_ZERO_KEY,
                VCFConstants.DBSNP_KEY,
                VCFConstants.SOMATIC_KEY);

        // TODO copy from TCGA spec..
        headerInfo.add(new VCFInfoHeaderLine("VT", 1, VCFHeaderLineType.String, "Variant type, can be SNP, INS or DEL"));


        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_ALLELE_DEPTHS,
                VCFConstants.GENOTYPE_PL_KEY);

        // cancer-specific
        // TODO: push to VCFConstants in GATK
        headerInfo.add(new VCFFormatHeaderLine("FA", VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Allele fraction of the alternate allele with regard to reference"));
        headerInfo.add(new VCFFormatHeaderLine("SS", 1, VCFHeaderLineType.Integer, "Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown"));
        headerInfo.add(new VCFFormatHeaderLine(VCFConstants.RMS_BASE_QUALITY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Average base quality for reads supporting alleles"));

        return headerInfo;
    }

    private VariantContext generateVC(CandidateMutation m) {
        GenomeLoc l = m.getLocation();
        List<Allele> alleles = Arrays.asList(Allele.create((byte) m.getRefAllele(), true), Allele.create((byte) m.getAltAllele()));
        List<Allele> tumorAlleles = Arrays.asList(Allele.create((byte) m.getRefAllele(), true), Allele.create((byte) m.getAltAllele()));
        List<Allele> normalAlleles = Arrays.asList(Allele.create((byte) m.getRefAllele(), true));

        GenotypeBuilder tumorGenotype =
                new GenotypeBuilder(m.getTumorSampleName(), tumorAlleles)
                        .AD(new int[]{m.getInitialTumorRefCounts(), m.getInitialTumorAltCounts()})
                        .attribute("FA", m.getTumorF())
                        .DP(m.getInitialTumorReadDepth());

        if (m.getInitialTumorAltCounts() > 0) {
            tumorGenotype.attribute(VCFConstants.RMS_BASE_QUALITY_KEY,  m.getInitialTumorAltQualitySum() / m.getInitialTumorAltCounts()); // TODO: is this TCGA compliant?
        }

        GenotypeBuilder normalGenotype =
                new GenotypeBuilder(m.getNormalSampleName(), normalAlleles)
                        .AD(new int[]{m.getInitialNormalRefCounts(), m.getInitialNormalAltCounts()})
                        .attribute("FA", m.getNormalF())
                        .DP(m.getInitialNormalReadDepth())
                        .attribute(VCFConstants.RMS_BASE_QUALITY_KEY, "." );    // TODO: is this TCGA-compliant?

        VariantContextBuilder vc =
                new VariantContextBuilder("", l.getContig(), l.getStart(), l.getStop(), alleles);

        vc.filter(m.isRejected()?"REJECT":"PASS");
        if(m.getDbsnpVC() != null) {
            vc.id(m.getDbsnpVC().getID());
            vc.attribute(VCFConstants.DBSNP_KEY, null);
        }
        if (!m.isRejected()) {
            vc.attribute(VCFConstants.SOMATIC_KEY, null);
            vc.attribute("VT", "SNP");
            tumorGenotype.attribute("SS", 2); // TODO: extract these TCGA specific attributes to a class
        }

        // add the genotype objects
        vc.genotypes(tumorGenotype.make(), normalGenotype.make());

        return vc.make();
    }

}
