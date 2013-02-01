package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMUtils;
import org.apache.commons.lang.math.NumberUtils;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.GenotypeBuilder;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;
import org.broadinstitute.sting.utils.variantcontext.writer.VariantContextWriter;

import java.io.PrintStream;
import java.util.*;

@PartitionBy(PartitionType.LOCUS)
@Reference(window=@Window(start=-1* MuTect.REFERENCE_HALF_WINDOW_LENGTH,stop= MuTect.REFERENCE_HALF_WINDOW_LENGTH))
@By(DataSource.REFERENCE)
public class MuTect extends LocusWalker<Integer, Integer> implements TreeReducible<Integer> {
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

    @Output(doc="VCF output of mutation candidates",shortName="vcf", fullName="vcf", required=false)
    protected VariantContextWriter vcf = null;

    /***************************************/
    // Reference Metadata inputs
    /***************************************/
    @Input(fullName="dbsnp", shortName = "dbsnp", doc="VCF file of DBSNP information", required=false)
    public List<RodBinding<VariantContext>> dbsnpRod = Collections.emptyList();

    @Input(fullName="cosmic", shortName = "cosmic", doc="VCF file of COSMIC sites", required=false)
    public List<RodBinding<VariantContext>> cosmicRod = Collections.emptyList();

    @Input(fullName="normal_panel", shortName = "normal_panel", doc="VCF file of sites observed in normal", required=false)
    public List<RodBinding<VariantContext>> normalPanelRod = Collections.emptyList();

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

    private boolean hasTumorBam = false;
    private boolean hasNormalBam = false;

    private double contaminantAlternateFraction;

    private TumorPowerCalculator tumorPowerCalculator;
    private NormalPowerCalculator normalNovelSitePowerCalculator;
    private NormalPowerCalculator normalDbSNPSitePowerCalculator;
    private TumorPowerCalculator strandArtifactPowerCalculator;

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
            MTAC.NORMAL_SAMPLE_NAME = "none";
        }

        this.contaminantAlternateFraction = Math.max(MTAC.MINIMUM_MUTATION_CELL_FRACTION, MTAC.FRACTION_CONTAMINATION);

        // coverage related initialization
        double powerConstantEps = Math.pow(10, -1 * (MTAC.POWER_CONSTANT_QSCORE/10));

        this.tumorPowerCalculator = new TumorPowerCalculator(powerConstantEps, MTAC.TUMOR_LOD_THRESHOLD, this.contaminantAlternateFraction);
        this.normalNovelSitePowerCalculator = new NormalPowerCalculator(powerConstantEps, MTAC.NORMAL_LOD_THRESHOLD);
        this.normalDbSNPSitePowerCalculator = new NormalPowerCalculator(powerConstantEps, MTAC.NORMAL_DBSNP_LOD_THRESHOLD);
        this.strandArtifactPowerCalculator = new TumorPowerCalculator(powerConstantEps, MTAC.STRAND_ARTIFACT_LOD_THRESHOLD, 0.0f);

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
        out.println("## muTect v1.0." + VERSION.split(" ")[1]);
        out.println(callStatsGenerator.generateHeader());

        // initialize the VCF output
        if (vcf != null) {
            Set<String> samples = new HashSet<String>();
            samples.add(MTAC.TUMOR_SAMPLE_NAME);
            samples.add(MTAC.NORMAL_SAMPLE_NAME);
            Set<VCFHeaderLine> headerInfo = VCFGenerator.getVCFHeaderInfo();
            vcf.writeHeader(new VCFHeader(headerInfo, samples));
        }

        lastTime = System.currentTimeMillis();
    }

    public static int MAX_INSERT_SIZE = 10000;
    private long totalReadsProcessed = 0;
    private long binReadsProcessed = 0;
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

            logger.info(String.format("[MUTECT] Processed %d reads in %d ms", totalReadsProcessed, elapsedTime));
        }

        // an optimization to speed things up when there is no coverage
        if ( !MTAC.FORCE_OUTPUT && numberOfReads == 0) { return -1; }

        final char upRef = Character.toUpperCase(ref.getBaseAsChar());
        String sequenceContext = SequenceUtils.createSequenceContext(ref, 3);

        try {


            // only process bases where the reference is [ACGT], because the FASTA for HG18 has N,M and R!
            if (upRef != 'A' && upRef != 'C' && upRef != 'G' && upRef != 'T') {
                return -1;
            }

            ArrayList<PileupElement> tumorPileupElements = new ArrayList<PileupElement>();
            ArrayList<PileupElement> normalPileupElements = new ArrayList<PileupElement>();

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
            double tumorPower = tumorPowerCalculator.cachingPowerCalculation(tumorBaseCount, MTAC.POWER_CONSTANT_AF);

            double normalPowerNoSNPPrior = normalNovelSitePowerCalculator.cachingPowerCalculation(normalBaseCount);
            double normalPowerWithSNPPrior = normalDbSNPSitePowerCalculator.cachingPowerCalculation(normalBaseCount);

            double normalPower = (germlineAtRisk)?normalPowerWithSNPPrior:normalPowerNoSNPPrior;

            double combinedPower = tumorPower*normalPower;
            if (!hasNormalBam) {
                combinedPower = tumorPower;
            }

            powerWriter.writeCoverage(rawContext, combinedPower);

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
                    logger.info(String.format("[MUTECT] Inspected %d potential candidates", candidatesInspected));
                }

                candidate.setInitialTumorAltCounts(tumorReadPile.qualitySums.getCounts(altAllele));
                candidate.setInitialTumorRefCounts(tumorReadPile.qualitySums.getCounts(upRef));
                candidate.setInitialTumorAltQualitySum(tumorReadPile.qualitySums.getQualitySum(altAllele));
                candidate.setInitialTumorRefQualitySum(tumorReadPile.qualitySums.getQualitySum(upRef));

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

                    contaminantLikelihoods.add(base, pe.getQual());
                }
                double[] refHetHom = LocusReadPile.extractRefHetHom(contaminantLikelihoods, upRef, altAllele);
                double contaminantLod = refHetHom[1] - refHetHom[0];
                candidate.setContaminantLod(contaminantLod);

                final QualitySums normQs = normalReadPile.qualitySums;


                VariableAllelicRatioGenotypeLikelihoods normalGl = normalReadPile.calculateLikelihoods(normalReadPile.qualityScoreFilteredPileup); // use MAPQ0 reads
                candidate.setInitialNormalBestGenotype(normalReadPile.getBestGenotype(normalGl));
                candidate.setInitialNormalLod(LocusReadPile.getRefVsAlt(normalGl, upRef, altAllele));


                double normalF = Math.max(LocusReadPile.estimateAlleleFraction(normalReadPile.qualityScoreFilteredPileup, upRef, altAllele), MTAC.MINIMUM_NORMAL_ALLELE_FRACTION);
                candidate.setNormalF(normalF);


                candidate.setInitialNormalAltQualitySum(normQs.getQualitySum(altAllele));
                candidate.setInitialNormalRefQualitySum(normQs.getQualitySum(upRef));

                candidate.setInitialNormalAltCounts(normQs.getCounts(altAllele));
                candidate.setInitialNormalRefCounts(normQs.getCounts(upRef));
                candidate.setInitialNormalReadDepth(normalReadPile.finalPileupReads.size());

                // TODO: parameterize filtering Mate-Rescued Reads (if someone wants to disable this)
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

                //TODO: shouldn't this be f2 in the lod calculation instead of the strand specific f values?
                //TODO: clean up use of forward/reverse vs positive/negative (prefer the latter since GATK uses it)
                ReadBackedPileup forwardPileup = filterReads(ref, tumorReadPile.finalPileupPositiveStrand, true).finalPileupPositiveStrand;
                double f2forward = LocusReadPile.estimateAlleleFraction(forwardPileup, upRef, altAllele);
                candidate.setTumorLodFStarForward(t2.calculateAltVsRefLOD(forwardPileup, (byte)altAllele, f2forward, 0.0));

                ReadBackedPileup reversePileup = filterReads(ref,tumorReadPile.finalPileupNegativeStrand, true).finalPileupNegativeStrand;
                double f2reverse = LocusReadPile.estimateAlleleFraction(reversePileup, upRef, altAllele);
                candidate.setTumorLodFStarReverse(t2.calculateAltVsRefLOD(reversePileup, (byte)altAllele, f2reverse, 0.0));

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

                candidate.setStrandContingencyTable(SequenceUtils.getStrandContingencyTable(refPile, mutantPile));

                // start with just the tumor pile
                candidate.setTumorAltForwardOffsetsInRead(SequenceUtils.getForwardOffsetsInRead(mutantPileup));
                candidate.setTumorAltReverseOffsetsInRead(SequenceUtils.getReverseOffsetsInRead(mutantPileup));

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

                if (MTAC.FORCE_ALLELES) {
                    out.println(callStatsGenerator.generateCallStats(candidate));
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
                        vcf.add(VCFGenerator.generateVC(m));
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

        if (candidate.isGermlineAtRisk() && candidate.getInitialNormalLod() < MTAC.NORMAL_DBSNP_LOD_THRESHOLD) {
            candidate.addRejectionReason("germline_risk");
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
    IndexedFastaSequenceFile refReader;

    private LocusReadPile filterReads(final ReferenceContext ref, final ReadBackedPileup pile, boolean filterMateRescueReads) {
        ArrayList<PileupElement> newPileupElements = new ArrayList<PileupElement>();

        for ( PileupElement p : pile ) {
            final GATKSAMRecord read = p.getRead();

            int mismatchQualitySum =
                    AlignmentUtils.mismatchesInRefWindow(p, ref, false, true);

            // do we have to many mismatches overall?
            if (mismatchQualitySum > this.MAX_READ_MISMATCH_QUALITY_SCORE_SUM) {
                continue;
            }

            // is this a heavily clipped read?
            if (SequenceUtils.isReadHeavilySoftClipped(read, MTAC.HEAVILY_CLIPPED_READ_FRACTION)) {
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


        return new LocusReadPile(newPileup, (char)ref.getBase(), 0, 0);
    }

    public enum ReadSource { Tumor, Normal }

    private ReadSource getReadSource(SAMRecord read) {
        // check if it's a tumor
        SAMReaderID id = getToolkit().getReaderIDForRead(read);
        if (tumorSAMReaderIDs.contains(id)) { return ReadSource.Tumor; }
        if (normalSAMReaderIDs.contains(id)) { return ReadSource.Normal; }

        // unexpected condition
        throw new RuntimeException("Unable to determine read source (tumor,normal) for read " + read.getReadName());
    }





}
