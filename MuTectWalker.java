package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import java.io.*;
import java.util.*;

import net.sf.picard.util.FormatUtil;
import net.sf.samtools.*;
import net.sf.samtools.util.StringUtil;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotype;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

@PartitionBy(PartitionType.INTERVAL)
@Reference(window=@Window(start=-1* MuTectWalker.REFERENCE_HALF_WINDOW_LENGTH,stop= MuTectWalker.REFERENCE_HALF_WINDOW_LENGTH))
@By(DataSource.REFERENCE)
public class MuTectWalker extends LocusWalker<Integer, Integer> implements TreeReducible<Integer> {
    private static final String TAB = "\t";
    public static final int REFERENCE_HALF_WINDOW_LENGTH = 150;
    public static final String BAM_TAG_TUMOR = "tumor";
    public static final String BAM_TAG_NORMAL = "normal";


    // DO NOT CHANGE THIS LINE!  It's the SVN revision number of the caller, which updates automatically!
    private static final String VERSION = "$Rev$";

    @Output(doc="Write output to here")
    PrintStream out;

    @Input(fullName="dbsnp", shortName = "dbsnp", doc="VCF file of DBSNP information", required=false)
    public RodBinding<VariantContext> dbsnpRod;

    @Input(fullName="cosmic", shortName = "cosmic", doc="VCF file of COSMIC sites", required=false)
    public RodBinding<VariantContext> cosmicRod;

    @Hidden
    @Input(fullName="normal_panel", shortName = "normal_panel", doc="VCF file of sites observed in normal", required=false)
    public RodBinding<VariantContext> normalPanelRod;

    @Argument(fullName = "noop", required = false, doc="used for debugging, basically exit as soon as we get the reads")
    public boolean NOOP = false;

    @Hidden
    @Argument(fullName = "enable_extended_output", required = false, doc="add many additional columns of statistics to the output file")
    public boolean ENABLE_EXTENDED_OUTPUT = false;

    @Hidden
    @Argument(fullName = "artifact_detection_mode", required = false, doc="used when running the caller on a normal (as if it were a tumor) to detect artifacts")
    public boolean ARTIFACT_DETECTION_MODE = false;

    @Argument(fullName = "tumor_sample_name", required = false, doc="name to use for tumor in output files")
    public String TUMOR_SAMPLE_NAME = null;

    @Argument(fullName = "bam_tumor_sample_name", required = false, doc="if the tumor bam contains multiple samples, only use read groups with SM equal to this value")
    public String BAM_TUMOR_SAMPLE_NAME = null;

    @Argument(fullName = "normal_sample_name", required = false, doc="name to use for normal in output files")
    public String NORMAL_SAMPLE_NAME = null;

    @Argument(fullName = "force_output", required = false, doc="force output for each site")
    public boolean FORCE_OUTPUT = false;

    @Argument(fullName = "force_alleles", required = false, doc="force output for all alleles at each site")
    public boolean FORCE_ALLELES = false;

    @Argument(fullName = "initial_tumor_lod", required = false, doc = "Initial LOD threshold for calling tumor variant")
    public float INITIAL_TUMOR_LOD_THRESHOLD = 6.3f;

    @Argument(fullName = "initial_tumor_fstar_lod", required = false, doc = "Initial F-star LOD threshold for calling tumor variant")
    public float INITIAL_TUMOR_FSTAR_LOD_THRESHOLD = 4.0f;

    @Argument(fullName = "tumor_lod", required = false, doc = "LOD threshold for calling tumor variant")
    public float TUMOR_LOD_THRESHOLD = 6.3f;

    @Argument(fullName = "fraction_contamination", required = false, doc = "estimate of fraction (0-1) of physical contamination with other unrelated samples")
    public float FRACTION_CONTAMINATION = 0.02f;

    @Argument(fullName = "minimum_mutation_cell_fraction", required = false,
            doc = "minimum fraction of cells which are presumed to have a mutation, used to handle non-clonality and contamination")
    public float MINIMUM_MUTATION_CELL_FRACTION = 0.00f;

    @Argument(fullName = "normal_lod", required = false, doc = "LOD threshold for calling normal non-variant")
    public float NORMAL_LOD_THRESHOLD = 2.3f;

    @Argument(fullName = "dbsnp_normal_lod", required = false, doc = "LOD threshold for calling normal non-variant at dbsnp sites")
    public float NORMAL_DBSNP_LOD_THRESHOLD = 5.3f;

    @Argument(fullName = "minimum_normal_allele_fraction", required = false, doc = "minimum allele fraction to be considered in normal, useful for normal sample contaminated with tumor")
    public float MINIMUM_NORMAL_ALLELE_FRACTION = 0.00f;

    @Argument(fullName = "tumor_f_pretest", required = false, doc = "for computational efficiency, reject sites with allelic fraction below this threshold")
    public float TUMOR_F_PRETEST = 0.005f;

    @Argument(fullName = "min_qscore", required = false, doc = "threshold for minimum base quality score")
    public int MIN_QSCORE = 5;

    public int MIN_QSUM_QSCORE = 13;

    @Argument(fullName = "gap_events_threshold", required = false, doc = "how many gapped events (ins/del) are allowed in proximity to this candidate")
    public int GAP_EVENTS_THRESHOLD = 3;

    @Argument(fullName = "heavily_clipped_read_fraction", required = false, doc = "if this fraction or more of the bases in a read are soft/hard clipped, do not use this read for mutation calling")
    public float HEAVILY_CLIPPED_READ_FRACTION = 0.30f;

    @Argument(fullName = "clipping_bias_pvalue_threshold", required = false, doc = "pvalue threshold for fishers exact test of clipping bias in mutant reads vs ref reads")
    public float CLIPPING_BIAS_PVALUE_THRESHOLD = 0.05f;

    @Argument(fullName = "min_mutant_sum_pretest", required = false, doc="pretest -- require sum of alt allele quality score to be >= this")
    public int MIN_MUTANT_SUM_PRETEST = 60;

    @Argument(fullName = "fraction_mapq0_threshold", required = false, doc = "threshold for determining if there is relatedness between the alt and ref allele read piles")
    public float FRACTION_MAPQ0_THRESHOLD = 0.5f;

    @Argument(fullName = "pir_median_threshold", required = false, doc="threshold for clustered read position artifact median")
    public double PIR_MEDIAN_THRESHOLD = 10;

    @Argument(fullName = "pir_mad_threshold", required = false, doc="threshold for clustered read position artifact MAD")
    public double PIR_MAD_THRESHOLD = 3;

    /** Parameters for ALT ALLELE IN NORMAL filter **/
    @Argument(fullName = "max_alt_alleles_in_normal_count", required = false, doc="threshold for maximum alternate allele counts in normal")
    public int MAX_ALT_ALLELES_IN_NORMAL_COUNT = 2;

    @Argument(fullName = "max_alt_alleles_in_normal_qscore_sum", required = false, doc="threshold for maximum alternate allele quality score sum in normal")
    public int MAX_ALT_ALLELES_IN_NORMAL_QSCORE_SUM = 20;

    @Argument(fullName = "max_alt_allele_in_normal_fraction", required = false, doc="threshold for maximum alternate allele fraction in normal")
    public double MAX_ALT_ALLELE_IN_NORMAL_FRACTION = 0.00;
    /******************************************/

    public boolean USE_MAPQ0_IN_NORMAL_QSCORE = true;

    /***************************************/
    // coverage related parameters
    /***************************************/
    @Argument(fullName="coverage_file", shortName="cov", doc="write out coverage in WIGGLE format to this file", required=false)
    public PrintStream COVERAGE_FILE = null;

    @Argument(fullName="coverage_20_q20_file", shortName="cov_q20", doc="write out 20x of Q20 coverage in WIGGLE format to this file", required=false)
    public PrintStream COVERAGE_20_Q20_FILE = null;

    @Argument(fullName="power_file", shortName="pow", doc="write out power in WIGGLE format to this file", required=false)
    public PrintStream POWER_FILE = null;

    @Argument(fullName="tumor_depth_file", shortName="tdf", doc="write out tumor read depth in WIGGLE format to this file", required=false)
    public PrintStream TUMOR_DEPTH_FILE = null;

    @Argument(fullName="normal_depth_file", shortName="ndf", doc="write out normal read depth in WIGGLE format to this file", required=false)
    public PrintStream NORMAL_DEPTH_FILE = null;

    @Argument(fullName="power_constant_qscore", doc="Phred scale quality score constant to use in power calculations", required=false)
    public int POWER_CONSTANT_QSCORE = 30;

    @Argument(fullName="absolute_copy_number_data", doc="Absolute Copy Number Data, as defined by Absolute, to use in power calculations", required=false)
    public File ABSOLUTE_COPY_NUMBER_DATA = null;

    @Argument(fullName="power_constant_af", doc="Allelic fraction constant to use in power calculations", required=false)
    public double POWER_CONSTANT_AF = 0.3f;

    public enum SequencingErrorModel {
        solid(5),
        illumina(1);

        private int priorBaseOffset;
        private SequencingErrorModel(int priorBaseOffset) {
            this.priorBaseOffset = priorBaseOffset;
        }

        public int getPriorBaseOffset() {
            return priorBaseOffset;
        }
    }

    @Hidden
    @Argument(fullName="sequencing_error_model", shortName="fP", required=false, doc="If provided, the error model will be forced to be the provided String. Valid options are illumina and solid (illumina works well as a generic model")
    public SequencingErrorModel SEQ_ERROR_MODEL = SequencingErrorModel.illumina;

    private static final FisherExact fisher = new FisherExact(5000);
    private final FormatUtil fmt = new FormatUtil();

    private boolean hasTumorBam = false;
    private boolean hasNormalBam = false;

    private double contaminantAlternateFraction;

    private double powerConstantEps;
    private TumorPowerCalculator tumorPowerCalculator;
    private NormalPowerCalculator normalNovelSitePowerCalculator;
    private NormalPowerCalculator normalDbSNPSitePowerCalculator;

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

    private int[] minimalHeaderIndicies;

    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    @Override
	public void initialize() {
        if (NOOP) { return; }

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
                    if (TUMOR_SAMPLE_NAME == null) {
                        try {
                            if (getToolkit().getReadsDataSource().getHeader(id).getReadGroups().size() == 0) {
                                throw new RuntimeException("No Read Groups found for Tumor BAM -- Read Groups are Required, or supply tumor_sample_name!");
                            }
                            TUMOR_SAMPLE_NAME = getToolkit().getReadsDataSource().getHeader(id).getReadGroups().get(0).getSample();
                        } catch (NullPointerException npe) {
                            TUMOR_SAMPLE_NAME = "tumor";
                        }
                    }
                } else if (BAM_TAG_NORMAL.equalsIgnoreCase(tag)) {
                    hasNormalBam = true;

                    // fill in the sample name if necessary
                    if (NORMAL_SAMPLE_NAME == null) {
                        try {
                            if (getToolkit().getReadsDataSource().getHeader(id).getReadGroups().size() == 0) {
                                throw new RuntimeException("No Read Groups found for Normal BAM -- Read Groups are Required, or supply normal_sample_name!");
                            }

                            NORMAL_SAMPLE_NAME = getToolkit().getReadsDataSource().getHeader(id).getReadGroups().get(0).getSample();
                        } catch (NullPointerException npe) {
                            NORMAL_SAMPLE_NAME = "normal";
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
            NORMAL_LOD_THRESHOLD = -1 * Float.MAX_VALUE;
            NORMAL_DBSNP_LOD_THRESHOLD = -1 * Float.MAX_VALUE;
        }

        this.contaminantAlternateFraction = Math.max(MINIMUM_MUTATION_CELL_FRACTION, FRACTION_CONTAMINATION);

        // coverage related initialization
        this.powerConstantEps = Math.pow(10, -1 * (POWER_CONSTANT_QSCORE/10));

        this.tumorPowerCalculator = new TumorPowerCalculator(this.powerConstantEps, TUMOR_LOD_THRESHOLD, this.contaminantAlternateFraction);
        this.normalNovelSitePowerCalculator = new NormalPowerCalculator(this.powerConstantEps, NORMAL_LOD_THRESHOLD);
        this.normalDbSNPSitePowerCalculator = new NormalPowerCalculator(this.powerConstantEps, NORMAL_DBSNP_LOD_THRESHOLD);

        stdCovWriter = new CoverageWiggleFileWriter(COVERAGE_FILE);
        q20CovWriter = new CoverageWiggleFileWriter(COVERAGE_20_Q20_FILE);
        powerWriter = new CoverageWiggleFileWriter(POWER_FILE);
        tumorDepthWriter = new CoverageWiggleFileWriter(TUMOR_DEPTH_FILE);
        normalDepthWriter = new CoverageWiggleFileWriter(NORMAL_DEPTH_FILE);

        // to force output, all we have to do is lower the initial tumor lod threshold to -infinity
        if (FORCE_OUTPUT) {
            this.INITIAL_TUMOR_LOD_THRESHOLD = -Float.MAX_VALUE;
        }

        // write out the call stats header
        out.println("## muTector v1.0." + VERSION.split(" ")[1]);
        String header;
        if (ENABLE_EXTENDED_OUTPUT) {
            header = StringUtil.join(TAB, COMPLETE_CALL_STATS_HEADER);
        } else {
            header = StringUtil.join(TAB, MINIMAL_CALL_STATS_HEADER);

            // initialize the indicies of the reduced headers from the full headers
            minimalHeaderIndicies = new int[MINIMAL_CALL_STATS_HEADER.length];
            for(int i=0; i<MINIMAL_CALL_STATS_HEADER.length; i++) {
                String column = MINIMAL_CALL_STATS_HEADER[i];

                for(int j=0; j<COMPLETE_CALL_STATS_HEADER.length; j++) {
                    if (COMPLETE_CALL_STATS_HEADER[j].equals(column)) {
                        minimalHeaderIndicies[i] = j;
                    }
                }
            }

        }
        out.println(header);

        lastTime = System.currentTimeMillis();
    }

    public static int MAX_INSERT_SIZE = 10000;
    private int totalReadsProcessed = 0;
    private int binReadsProcessed = 0;
    private long lastTime;
    private int candidatesInspected = 0;

    @Override
	public Integer map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext rawContext) {
        if (NOOP) return 0;
        
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
        if ( !FORCE_OUTPUT && numberOfReads == 0) { return -1; }

        final char upRef = Character.toUpperCase(ref.getBaseAsChar());
        try {


            // only process bases where the reference is [ACGT], because the FASTA for HG18 has N,M and R!
            if (upRef != 'A' && upRef != 'C' && upRef != 'G' && upRef != 'T') {
//                final String msg = "FAILED non-ACGT reference.  Reference was " + upRef;
//                writeCallStats(CallStatsVerbosity.MINOR, rawContext, msg);
                return -1;
            }

            ArrayList<PileupElement> normalPileupElements = new ArrayList<PileupElement>();
            ArrayList<PileupElement> tumorPileupElements = new ArrayList<PileupElement>();

            int totalPairs = 0;
            int improperPairs = 0;
            List<PileupElement> primaryPileup = new ArrayList<PileupElement>(pileup.depthOfCoverage());
            List<PileupElement> secondaryPileup = new ArrayList<PileupElement>(pileup.depthOfCoverage());
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

            if (BAM_TUMOR_SAMPLE_NAME != null) {
                tumorPileup = tumorPileup.getPileupForSample(BAM_TUMOR_SAMPLE_NAME);
            }

            // TODO: do this filtering much earlier so that in the simulation we are drawing from reads that will pass this filter!
            final LocusReadPile tumorReadPile = new LocusReadPile(tumorPileup, upRef, MIN_QSCORE, MIN_QSUM_QSCORE, false, ARTIFACT_DETECTION_MODE);
            final LocusReadPile normalReadPile = new LocusReadPile(normalPileup, upRef ,MIN_QSCORE, 0, this.USE_MAPQ0_IN_NORMAL_QSCORE, true);


            Collection<VariantContext> panelOfNormalsVC = tracker.getValues(normalPanelRod, rawContext.getLocation());
            Collection<VariantContext> cosmicVC = tracker.getValues(cosmicRod, rawContext.getLocation());
            Collection<VariantContext> dbsnpVC = tracker.getValues(dbsnpRod, rawContext.getLocation());

            // remove the effect of cosmic from dbSNP
            boolean knownDbSnpSite = (!dbsnpVC.isEmpty() && cosmicVC.isEmpty());

            // compute coverage flags
            int tumorCoveredDepthThreshold = 14;
            int normalCoveredDepthThreshold = (knownDbSnpSite)?19:8;
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
            double combinedPower;
            if (ABSOLUTE_COPY_NUMBER_DATA == null) {
                tumorPower = tumorPowerCalculator.cachingPowerCalculation(tumorBaseCount, POWER_CONSTANT_AF);

                NormalPowerCalculator npc = (knownDbSnpSite)?normalDbSNPSitePowerCalculator:normalNovelSitePowerCalculator;
                normalPower = npc.cachingPowerCalculation(normalBaseCount);

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
                if (!FORCE_OUTPUT && tumorReadPile.qualitySums.getCounts(altAllele) == 0) { continue; }

                CandidateMutation candidate = new CandidateMutation(rawContext.getLocation(), upRef);
                candidate.setTumorSampleName(TUMOR_SAMPLE_NAME);
                candidate.setNormalSampleName(NORMAL_SAMPLE_NAME);
                candidate.setCovered(isBaseCovered);
                candidate.setPower(combinedPower);
                candidate.setTumorPower(tumorPower);
                candidate.setNormalPower(normalPower);
                candidate.setTumorQ20Count(tumorQ20BaseCount);
                candidate.setNormalQ20Count(normalQ20BaseCount);
                candidate.setInitialTumorNonRefQualitySum(tumorReadPile.qualitySums.getOtherQualities(upRef));
                candidate.setAltAllele(altAllele);
                candidate.setTotalPairs(totalPairs);
                candidate.setImproperPairs(improperPairs);
                candidate.setMapQ0Reads(mapQ0Reads);
                candidate.setEstimatedFractionContamination(FRACTION_CONTAMINATION);
                candidate.setPanelOfNormalsVC(panelOfNormalsVC.isEmpty()?null:panelOfNormalsVC.iterator().next());
                candidate.setCosmicSite(!cosmicVC.isEmpty());
                candidate.setDbsnpSite(knownDbSnpSite);

                candidate.setTumorF(tumorReadPile.estimateAlleleFraction(upRef, altAllele));
                if (!FORCE_OUTPUT && candidate.getTumorF() < TUMOR_F_PRETEST) {
                    continue;
                }

                if (++candidatesInspected % 1000 == 0) {
                    logger.info(String.format("[MUTECTOR] Inspected %d potential candidates", candidatesInspected));
                }
                
                candidate.setInitialTumorAltCounts(tumorReadPile.qualitySums.getCounts(altAllele));
                candidate.setInitialTumorRefCounts(tumorReadPile.qualitySums.getCounts(upRef));
                candidate.setInitialTumorAltQualitySum(tumorReadPile.qualitySums.getQualitySum(altAllele));
                candidate.setInitialTumorRefQualitySum(tumorReadPile.qualitySums.getQualitySum(upRef));

                VariableAllelicRatioGenotypeLikelihoods tumorGl = tumorReadPile.calculateLikelihoods(tumorReadPile.finalPileup);
                candidate.setInitialTumorLod(tumorReadPile.getAltVsRef(tumorGl, upRef, altAllele));
                candidate.setInitialTumorReadDepth(tumorReadPile.finalPileupReads.size());

                VariableAllelicRatioGenotypeLikelihoods tumorFStarGl = tumorReadPile.calculateLikelihoods(candidate.getTumorF(), tumorReadPile.finalPileup);
                candidate.setTumorLodFStar(tumorReadPile.getHetVsRef(tumorFStarGl, upRef, altAllele));

                candidate.setTumorInsertionCount(tumorReadPile.getInsertionsCount());
                candidate.setTumorDeletionCount(tumorReadPile.getDeletionsCount());

                if (candidate.getInitialTumorLod() >= this.INITIAL_TUMOR_LOD_THRESHOLD ||
                    candidate.getTumorLodFStar() >= this. INITIAL_TUMOR_FSTAR_LOD_THRESHOLD ) {

                        // if only Java had "unless" and I didn't hate deMorgan's so much...
                } else {
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
                    contaminantLikelihoods.add(base, pe.getQual());
                }
                double[] refHetHom = LocusReadPile.extractRefHetHom(contaminantLikelihoods, upRef, altAllele);
                double contaminantLod = refHetHom[1] - refHetHom[0];
                candidate.setContaminantLod(contaminantLod);

//                // perform rank-sum test on distribution of quality scores
//                // for tumor reference allele and alternate allele
//                int[] refQuals = convertToSortedArray(tumorReadPile.getQualityScores(upRef));
//                int[] altQuals = convertToSortedArray(tumorReadPile.getQualityScores(altAllele));
//
//                // only works if we have both ref and alt alleles (e.g. not hnref case!)
//                if (refQuals.length > 0 && altQuals.length > 0) {
//                    candidate.setTumorQualityRankSumTest(getRankSumTest().test(refQuals, altQuals));
//                }
//


                // (ii) the quality score sum for the mutant base in the normal must be < 50 and the
                //      LOD score for ref:ref vs mutant:ref + mutant:mutant must be at least 2.3.
                final QualitySums normQs = normalReadPile.qualitySums;

                VariableAllelicRatioGenotypeLikelihoods normalGl = normalReadPile.calculateLikelihoods(normalReadPile.qualityScoreFilteredPileup); // use MAPQ0 reads
                candidate.setInitialNormalBestGenotype(normalReadPile.getBestGenotype(normalGl));
                candidate.setInitialNormalLod(LocusReadPile.getRefVsAlt(normalGl, upRef, altAllele));

                double normalF = Math.max(normalReadPile.estimateAlleleFraction(normalReadPile.qualityScoreFilteredPileup, upRef, altAllele), MINIMUM_NORMAL_ALLELE_FRACTION);
                candidate.setNormalF(normalF);
                VariableAllelicRatioGenotypeLikelihoods normalFStarGl = normalReadPile.calculateLikelihoods(normalF, normalReadPile.qualityScoreFilteredPileup);
                candidate.setNormalLodFStar(normalReadPile.getRefVsHet(normalFStarGl, upRef, altAllele));

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
                    priorBasePositiveDirection = (char) ref.getBases()[(int)candidate.getLocation().getStart() - SEQ_ERROR_MODEL.getPriorBaseOffset() - (int)ref.getWindow().getStart()];
                }

                if (containsPosition(ref.getWindow(), candidate.getLocation().getStart() + 1)) {
                    priorBaseNegativeDirection = (char) ref.getBases()[(int)candidate.getLocation().getStart() + SEQ_ERROR_MODEL.getPriorBaseOffset() - (int)ref.getWindow().getStart()];
                }

                candidate.setPriorBasePositiveDirection(priorBasePositiveDirection);
                candidate.setPriorBaseNegativeDirection(priorBaseNegativeDirection);

                final LocusReadPile t2 = filterReads(ref, tumorReadPile, altAllele, refGATKString, refStart);

                // if there are no reads remaining, abandon this theory
                if ( !FORCE_OUTPUT && t2.finalPileupReads.size() == 0) { continue; }

                candidate.setInitialTumorAltCounts(t2.qualitySums.getCounts(altAllele));
                candidate.setInitialTumorRefCounts(t2.qualitySums.getCounts(upRef));
                candidate.setInitialTumorAltQualitySum(t2.qualitySums.getQualitySum(altAllele));
                candidate.setInitialTumorRefQualitySum(t2.qualitySums.getQualitySum(upRef));

                VariableAllelicRatioGenotypeLikelihoods t2Gl = t2.calculateLikelihoods(t2.finalPileup);
                candidate.setInitialTumorLod(t2.getAltVsRef(t2Gl, upRef, altAllele));
                candidate.setInitialTumorReadDepth(t2.finalPileupReads.size());

                double f2 = t2.estimateAlleleFraction(upRef, altAllele);
                VariableAllelicRatioGenotypeLikelihoods tumorFStarGl2 = t2.calculateLikelihoods(f2, t2.finalPileup);
                candidate.setTumorLodFStar(tumorReadPile.getHetVsRef(tumorFStarGl2, upRef, altAllele));
                candidate.setTumorF(f2);

                ReadBackedPileup forwardPileup = t2.finalPileup.getPositiveStrandPileup();
                double f2forward = LocusReadPile.estimateAlleleFraction(forwardPileup, upRef, altAllele);
                VariableAllelicRatioGenotypeLikelihoods tumorFStarGl2Forward = t2.calculateLikelihoods(f2forward, forwardPileup);
                candidate.setTumorLodFStarForward(tumorReadPile.getHetVsRef(tumorFStarGl2Forward, upRef, altAllele));

                ReadBackedPileup reversePileup = t2.finalPileup.getNegativeStrandPileup();
                double f2reverse = LocusReadPile.estimateAlleleFraction(reversePileup, upRef, altAllele);
                VariableAllelicRatioGenotypeLikelihoods tumorFStarGl2Reverse = t2.calculateLikelihoods(f2reverse, reversePileup);
                candidate.setTumorLodFStarReverse(tumorReadPile.getHetVsRef(tumorFStarGl2Reverse, upRef, altAllele));

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

                // FIXME: shouldn't this be refAllele here?
                final LocusReadPile mutantPile = new LocusReadPile(mutantPileup, altAllele, 0, 0);
                final LocusReadPile refPile =  new LocusReadPile(referencePileup, altAllele, 0, 0);


                candidate.setPerfectStrandBias(calculatePerfectStrandBias(mutantPile));
                candidate.setStrandBias(calculateStrandBias(refPile, mutantPile));

                // start with just the tumor pile
                candidate.setTumorAltForwardOffsetsInRead(getForwardOffsetsInRead(mutantPileup));
                candidate.setTumorAltReverseOffsetsInRead(getReverseOffsetsInRead(mutantPileup));


//                candidate.setClassicSkewScoresAndOffsets(
//                    performClassicMisalignmentTest(refPile, mutantPile, rawContext.getContig(), rawContext.getPosition(), refGATKString, refStart, MIN_QSCORE)
//                );
//
//                candidate.setFisherSkewScoresAndOffsets(
//                    performFisherMisalignmentTest(refPile, mutantPile, rawContext.getContig(), rawContext.getPosition(), refGATKString, refStart, MIN_QSCORE)
//                );

                //candidate.setClippingBias(calculateClippingBias(refPile, mutantPile));


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

                if (FORCE_ALLELES) {
                    out.println(generateCallStats(candidate));
                } else {
                    messageByTumorLod.put(candidate.getInitialTumorLod(), generateCallStats(candidate));
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

    private boolean containsPosition(GenomeLoc window, int position) {
        return (window.getStart() <= position && position <= window.getStop());
    }

    // FIXME: there must be a better way to do this...
    private int[] convertToSortedArray(List inputList) {
        Collections.sort(inputList);
        int[] output = new int[inputList.size()];
        for(int i=0; i< inputList.size(); i++) {
            Object o = inputList.get(i);
            if (o instanceof Integer) {
                output[i] = (Integer) o;
            } else if (o instanceof Byte) {
                output[i] = ((Byte) o).intValue();
            } else {
                throw new RuntimeException("Unknown value for RST");
            }
        }
        return output;
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
            if (isReadHeavilyClipped(rec)) { b++;} else { a++; }
        }
        for (SAMRecord rec : mutantPile.finalPileupReads) {
            if (isReadHeavilyClipped(rec)) { d++;} else { c++; }
        }

        final double p = fisher.getTwoTailedP(a,b,c,d);
        int n = c+d;
        final double pmaxpos = fisher.getTwoTailedP(a,b,n,0);
        final double pmaxneg = fisher.getTwoTailedP(a,b,0,n);

        return new FisherData(a,b,c,d,p, pmaxpos, pmaxneg);
    }

    private boolean isReadHeavilyClipped(SAMRecord rec) {
        int total = 0;
        int clipped = 0;
        for(CigarElement ce : rec.getCigar().getCigarElements()) {
            total += ce.getLength();
            if (ce.getOperator() == CigarOperator.HARD_CLIP || ce.getOperator() == CigarOperator.SOFT_CLIP) {
                clipped += ce.getLength();
            }
        }

        return ((float) clipped / (float)total >= HEAVILY_CLIPPED_READ_FRACTION);
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
            // FIXME: maybe we should be doing start-site distribution, or a clipping aware offset?
                positions.add(
                        Math.abs((int)(p.getLocation().getStart() - (useForwardOffsets?pe.getRead().getAlignmentStart():pe.getRead().getAlignmentEnd())))
                );
        }

        return positions;
    }

    private static final double PERFECT_STRAND_BIAS_THRESHOLD = .05d;
    private static final double STRAND_BIAS_THRESHOLD = .05d;

    private static final int QUALITY_RST_MEDIAN_SHIFT_THRESHOLD = 4;
    private static final double QUALITY_RST_PVALUE_THRESHOLD = .005;

    private void performRejection(CandidateMutation candidate) {
        if (candidate.getTumorLodFStar() < TUMOR_LOD_THRESHOLD) {
            candidate.addRejectionReason("fstar_tumor_lod");
        }

        if (ARTIFACT_DETECTION_MODE) {
            return;
        }

        if (candidate.getInitialTumorAltQualitySum() < MIN_MUTANT_SUM_PRETEST) {
            candidate.addRejectionReason("min_mutant_sum_pretest");
        }

        // if the best theory for the normal is A het (not necessarily this het)
        // just move on.  We're not attempting to call LOH with this tool
        if (candidate.getInitialNormalBestGenotype() != DiploidGenotype.AA &&
            candidate.getInitialNormalBestGenotype() != DiploidGenotype.CC &&
            candidate.getInitialNormalBestGenotype() != DiploidGenotype.GG &&
            candidate.getInitialNormalBestGenotype() != DiploidGenotype.TT) {
            candidate.addRejectionReason("het_normal");
        }

        if (candidate.getTumorInsertionCount() >= GAP_EVENTS_THRESHOLD ||
            candidate.getTumorDeletionCount()  >= GAP_EVENTS_THRESHOLD) {
            candidate.addRejectionReason("nearby_gap_events");
        }

        if (FRACTION_CONTAMINATION+MINIMUM_MUTATION_CELL_FRACTION > 0 && candidate.getTumorLodFStar() <= TUMOR_LOD_THRESHOLD + Math.max(0, candidate.getContaminantLod())) {
            candidate.addRejectionReason("possible_contamination");
        }

        if (candidate.isDbsnpSite() && candidate.getInitialNormalLod() < NORMAL_DBSNP_LOD_THRESHOLD) {
            candidate.addRejectionReason("DBSNP Site");
        }

        if (candidate.getInitialNormalLod() < NORMAL_LOD_THRESHOLD) {
            candidate.addRejectionReason("normal_lod");
        }

        if ( (candidate.getInitialNormalAltCounts() >= MAX_ALT_ALLELES_IN_NORMAL_COUNT && candidate.getInitialNormalAltQualitySum() > MAX_ALT_ALLELES_IN_NORMAL_QSCORE_SUM && candidate.getNormalF() > MAX_ALT_ALLELE_IN_NORMAL_FRACTION)) {
            candidate.addRejectionReason("alt allele in normal");
        }

        if ( (candidate.getTumorForwardOffsetsInReadMedian() != null && candidate.getTumorForwardOffsetsInReadMedian() <= PIR_MEDIAN_THRESHOLD && candidate.getTumorForwardOffsetsInReadMad() != null && candidate.getTumorForwardOffsetsInReadMad() <= PIR_MAD_THRESHOLD) ||
              candidate.getTumorReverseOffsetsInReadMedian() != null && candidate.getTumorReverseOffsetsInReadMedian() <= PIR_MEDIAN_THRESHOLD && candidate.getTumorReverseOffsetsInReadMad() != null && candidate.getTumorReverseOffsetsInReadMad() <= PIR_MAD_THRESHOLD ) {
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
        if (candidate.isPositiveDirectionAtRisk() && candidate.isPositiveDirectionPowered() && candidate.getTumorLodFStarReverse() < 2.0) {
            candidate.addRejectionReason("positive_strand_artifact");
        }

        if (candidate.isNegativeDirectionAtRisk() && candidate.isNegativeDirectionPowered() && candidate.getTumorLodFStarForward() < 2.0) {
            candidate.addRejectionReason("negative_strand_artifact");
        }

        if (candidate.getTotalPairs() > 0 && ((float)candidate.getMapQ0Reads() / (float)candidate.getTotalPairs()) >= FRACTION_MAPQ0_THRESHOLD) {
            candidate.addRejectionReason("poor_mapping_region_mapq0");
        }


        if (candidate.isSeenInPanelOfNormals()) {
            candidate.addRejectionReason("seen_in_panel_of_normals");
        }
    }

    private Integer getKeyForSmallestValue(Map<Integer, Double> map) {
        Integer key = null;
        Double smallest = null;
        for (Map.Entry<Integer, Double> e : map.entrySet()) {
            if (smallest == null || e.getValue() < smallest) {
                key = e.getKey();
                smallest = e.getValue();
            }
        }
        return key;
    }

    private Integer getKeyForLargestValue(Map<Integer, Double> map) {
        Integer key = null;
        Double largest = null;
        for (Map.Entry<Integer, Double> e : map.entrySet()) {
            if (largest == null || e.getValue() > largest) {
                key = e.getKey();
                largest = e.getValue();
            }
        }
        return key;
    }

    private static final String[] COMPLETE_CALL_STATS_HEADER =
            new String[]{
                "contig","position","ref_allele","alt_allele","tumor_name","normal_name","score","dbsnp_site",
                "covered", "power", "tumor_power", "normal_power",
                "total_pairs","improper_pairs","map_Q0_reads",
                "init_t_lod","t_lod_fstar","t_lod_fstar_forward", "t_lod_fstar_reverse", "tumor_f","contaminant_fraction","contaminant_lod","minimum_tumor_f", "t_q20_count", "t_ref_count","t_alt_count","t_ref_sum","t_alt_sum","t_ins_count","t_del_count",
                "normal_best_gt","init_n_lod","n_lod_fstar","normal_f","n_q20_count", "n_ref_count","n_alt_count","n_ref_sum","n_alt_sum",
                "at_risk_positive_direction_artifact", "at_risk_negative_direction_artifact", "powered_positive_direction_artifact", "powered_negative_direction_artifact",
                "perfect_strand_bias","strand_bias_counts","strand_bias",
                "classic_max_skew_lod", "classic_max_skew_lod_offset", "fisher_min_skew_pvalue", "fisher_min_skew_pvalue_offset",
                "tumor_qsrst_ms", "tumor_qsrst_pval", "tumor_rprst_ms", "tumor_rprst_pval",
                "tumor_alt_fpir_median", "tumor_alt_fpir_mad","tumor_alt_rpir_median","tumor_alt_rpir_mad","alt_fpir","alt_rpir",
                "powered_filters", "observed_in_normals_count", "failure_reasons","judgement"
            };

    private static final String[] MINIMAL_CALL_STATS_HEADER =
            new String[]{
                "contig","position","ref_allele","alt_allele","tumor_name","normal_name","score","dbsnp_site",
                "covered", "power", "tumor_power", "normal_power",
                "total_pairs","improper_pairs","map_Q0_reads",
                "t_lod_fstar","tumor_f","contaminant_fraction","contaminant_lod",
                "t_ref_count","t_alt_count","t_ref_sum","t_alt_sum","t_ins_count","t_del_count",
                "normal_best_gt","init_n_lod", "n_ref_count","n_alt_count","n_ref_sum","n_alt_sum",
                "judgement"
            };

    private String generateCallStats(CandidateMutation candidate) {
        RankSumTest.Result qrst = candidate.getTumorQualityRankSumTest();
        RankSumTest.Result prst = candidate.getTumorReadPositionRankSumTest();

        Double classicSkewScore = null;
        Integer classicSkewOffset = null;
        Map<Integer, Double> classicSkewInfo = candidate.getClassicSkewScoresAndOffsets();
        if (classicSkewInfo != null && classicSkewInfo.size() > 0) {
            classicSkewOffset = getKeyForLargestValue(classicSkewInfo);
            classicSkewScore = classicSkewInfo.get(classicSkewOffset);
            if (classicSkewScore == Double.POSITIVE_INFINITY) {
                classicSkewScore = 999999d;
            } else if (classicSkewScore == Double.NEGATIVE_INFINITY) {
                classicSkewScore = -999999d;
            }

        }

        Double fisherSkewScore = null;
        Integer fisherSkewOffset = null;
        Map<Integer, Double> fisherSkewInfo = candidate.getFisherSkewScoresAndOffsets();
        if (fisherSkewInfo != null && fisherSkewInfo.size() > 0) {
            fisherSkewOffset = getKeyForSmallestValue(fisherSkewInfo);
            fisherSkewScore = fisherSkewInfo.get(fisherSkewOffset);
        }

        // further classify KEEP to indicate KEEP-CLASSIC for classic LOD
        String keepString = "REJECT";
        if (!candidate.isRejected()) {
            keepString = "KEEP";
        }

        String siteInfo = "NOVEL";
        if (candidate.isDbsnpSite()) {
            siteInfo = "DBSNP";
        }
        if (candidate.isCosmicSite()) {
            siteInfo = "COSMIC";
        }

        String[] msg = new String[] {
                        candidate.getLocation().getContig(),
                        format(candidate.getLocation().getStart()),
                        ""+candidate.getRefAllele(),
                        ""+candidate.getAltAllele(),
                        candidate.getTumorSampleName(),
                        candidate.getNormalSampleName(),
                        format(candidate.getScore()),
                        siteInfo,
                        (candidate.isCovered()?"COVERED":"UNCOVERED"),
                        format(candidate.getPower()),
                        format(candidate.getTumorPower()),
                        format(candidate.getNormalPower()),
                        format(candidate.getTotalPairs()),
                        format(candidate.getImproperPairs()),
                        format(candidate.getMapQ0Reads()),
                        format(candidate.getInitialTumorLod()),
                        format(candidate.getTumorLodFStar()),
                        format(candidate.getTumorLodFStarForward()),
                        format(candidate.getTumorLodFStarReverse()),
                        format(candidate.getTumorF()),
                        format(FRACTION_CONTAMINATION),
                        format(candidate.getContaminantLod()),
                        format("n/a"),
                        format(candidate.getTumorQ20Count()),
                        format(candidate.getInitialTumorRefCounts()),
                        format(candidate.getInitialTumorAltCounts()),
                        format(candidate.getInitialTumorRefQualitySum()),
                        format(candidate.getInitialTumorAltQualitySum()),
                        format(candidate.getTumorInsertionCount()),
                        format(candidate.getTumorDeletionCount()),
                        format(candidate.getInitialNormalBestGenotype().toString()),
                        format(candidate.getInitialNormalLod()),
                        format(candidate.getNormalLodFStar()),
                        format(candidate.getNormalF()),
                        format(candidate.getNormalQ20Count()),
                        format(candidate.getInitialNormalRefCounts()),
                        format(candidate.getInitialNormalAltCounts()),
                        format(candidate.getInitialNormalRefQualitySum()),
                        format(candidate.getInitialNormalAltQualitySum()),
                        format(candidate.isPositiveDirectionAtRisk()?1:0),
                        format(candidate.isNegativeDirectionAtRisk()?1:0),
                        format(candidate.isPositiveDirectionPowered()?1:0),
                        format(candidate.isNegativeDirectionPowered()?1:0),
                        format(candidate.getPerfectStrandBias().getP()),
                        format(candidate.getStrandBias().dataToString()),
                        format(candidate.getStrandBias().getP()),
                        format(classicSkewScore),
                        format(classicSkewOffset),
                        format(fisherSkewScore),
                        format(fisherSkewOffset),
                        qrst==null?"n/a":format(candidate.getTumorQualityRankSumTest().getMedianShift()),
                        qrst==null?"n/a":format(candidate.getTumorQualityRankSumTest().getP()),
                        prst==null?"n/a":format(candidate.getTumorReadPositionRankSumTest().getMedianShift()),
                        prst==null?"n/a":format(candidate.getTumorReadPositionRankSumTest().getP()),
                        candidate.getTumorForwardOffsetsInReadMedian()==null?"n/a":format(candidate.getTumorForwardOffsetsInReadMedian()),
                        candidate.getTumorForwardOffsetsInReadMad()==null?"n/a":format(candidate.getTumorForwardOffsetsInReadMad()),
                        candidate.getTumorReverseOffsetsInReadMedian()==null?"n/a":format(candidate.getTumorReverseOffsetsInReadMedian()),
                        candidate.getTumorReverseOffsetsInReadMad()==null?"n/a":format(candidate.getTumorReverseOffsetsInReadMad()),
                        (candidate.getTumorAltForwardOffsetsInRead().size() > 500 ? "too_many":StringUtil.join(",", toString(candidate.getTumorAltForwardOffsetsInRead()))),
                        (candidate.getTumorAltReverseOffsetsInRead().size() > 500 ? "too_many":StringUtil.join(",", toString(candidate.getTumorAltReverseOffsetsInRead()))),
                        StringUtil.join(",", candidate.getPoweredFilters().toArray(new String[]{})),
                        format(candidate.getCountOfNormalsObservedIn()),
                        StringUtil.join(",", candidate.getRejectionReasons().toArray(new String[]{})),
                        keepString
            };

        if (ENABLE_EXTENDED_OUTPUT) {
            return StringUtil.join(TAB, msg);
        } else {
            List<String> output = new ArrayList<String>();
            for(int index : minimalHeaderIndicies) {
                output.add(msg[index]);
            }
            return StringUtil.join(TAB, output.toArray(new String[]{}));
        }
    }

    private String format(String s) { return s; }
    private String format(Integer i) { return fmt.format(i); }
    private String format(Float f) { return format((double)f);}
    private String format(Double d) {
        if (d == null) { return "n/a"; }

        String s = fmt.format(d);
        return ("-0".equals(s))?"0":s;
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

    private String[] toString(List<Integer> ints) {
        String[] out = new String[ints.size()];
        for(int i=0; i<ints.size(); i++) {
            out[i] = ints.get(i).toString();
        }
        return out;
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

    @Override
    public void onTraversalDone(final Integer result) {
//        quietClose(this.mafWriter);
//        quietClose(this.coverageWriter);
    }

    // TODO: commented out since it has to be upgraded to handle fragment based calling, but if we no longer use this test then yay! :)

//    private Map<Integer, Double> performClassicMisalignmentTest(final LocusReadPile refPile, final LocusReadPile mutantPile, final String contig, final long position, final String reference, final long leftmostIndex, final int minQScore) {
//        final Map<Integer, Double> skewLodOffsets = new TreeMap<Integer, Double>();
//
//        //fixme: calculate this range properly
//        int MAX_OFFSET_DISTANCE = 60;
//
//        for(int offset=-1 * MAX_OFFSET_DISTANCE; offset<MAX_OFFSET_DISTANCE; offset++) {
//            // what is the reference at this position?
//            long refOffset = position - leftmostIndex + offset;
//
//            // cap it at the ends of the provided reference sequence also!
//            if (refOffset < 0 || refOffset >= reference.length()) { continue; }
//
//            // if this is a dbsnp site... exclude it also
//            if (isDbsnpSite(contig, (int)position + offset)) { continue; }
//
//            // allow for doubletons
//            if (offset >= -1 && offset <= 1 ) { continue; }
//
//            final int SKEW_QSCORE_THRESHOLD = 10;
//            final double[] mutantNormProbs = mutantPile.getNormalizedProbs(offset, SKEW_QSCORE_THRESHOLD);
//            final double[] otherNormProbs = refPile.getNormalizedProbs(offset, SKEW_QSCORE_THRESHOLD);
//
//            double J = 0;
//            for(int i=0; i<10; i++) {
//                J += mutantNormProbs[i] * otherNormProbs[i];
//            }
//            final double skewLod = Math.log10( (1-J) / J);
//
//
//            // FIXME!: get rid of bases with quality score < min Qscore (add parameter to getLocusBases)
//            List<Byte> mutantAlleles = mutantPile.getLocusBases(offset);
//            List<Byte> otherAlleles = refPile.getLocusBases(offset);
//
//            final int mutantReadCounts = mutantAlleles.size();
//            final int otherReadCounts = otherAlleles.size();
//
//            if (mutantReadCounts >= 2 && otherReadCounts >= 2) {
//                skewLodOffsets.put(offset, skewLod);
//            }
//        }
//
//        return skewLodOffsets;
//
//    }

//    private Map<Integer, Double> performFisherMisalignmentTest(final LocusReadPile refPile,
//                                                               final LocusReadPile mutantPile,
//                                                               final String contig,
//                                                               final long position,
//                                                               final String reference,
//                                                               final long leftmostIndex,
//                                                               final int minQScore) {
//        final Map<Integer, Double> skewLodOffsets = new TreeMap<Integer, Double>();
//
//        //fixme: calculate this range properly
//        int MAX_OFFSET_DISTANCE = 60;
//
//        for(int offset=-1 * MAX_OFFSET_DISTANCE; offset<MAX_OFFSET_DISTANCE; offset++) {
//
//            // what is the reference at this position?
//            long offsetIntoReferenceString = position - leftmostIndex + offset;
//
//            // cap it at the ends of the provided reference sequence also!
//            if (offsetIntoReferenceString < 0 || offsetIntoReferenceString >= reference.length()) { continue; }
//
//            // if this is a dbsnp site... exclude it also
//            if (isDbsnpSite(contig, (int) position + offset)) { continue; }
//
//            // allow for doubletons
//            if (offset >= -1 && offset <= 1 ) { continue; }
//
//
//            // FIXME!: get rid of bases with quality score < min Qscore (add parameter to getLocusBases)
//            List<Byte> mutantAlleles = mutantPile.getLocusBases(offset);
//            List<Byte> otherAlleles = refPile.getLocusBases(offset);
//
//            // determine the two most common alleles at the offset
//            TreeMap<Byte, Integer> counts = new TreeMap<Byte, Integer>();
//            counts.put((byte)'A',0);
//            counts.put((byte)'C',0);
//            counts.put((byte)'G',0);
//            counts.put((byte)'T',0);
//            for (byte base : otherAlleles) {
//                counts.put(base, (counts.containsKey(base)?counts.get(base):0) + 1);
//            }
//            for (byte base : mutantAlleles) {
//                counts.put(base, (counts.containsKey(base)?counts.get(base):0) + 1);
//            }
//
//            LinkedHashMap sortedCounts = (LinkedHashMap<Byte, Integer>) sortByDescendingValue(counts);
//            Iterator<Byte> it = sortedCounts.keySet().iterator();
//            byte altAllele1 = it.next();
//            byte altAllele2 = it.next();
//
//
//
//            // Construct a 2x2 contingency table of
//            //               OFFSET-ALT1   OFFSET-ALT2
//            //    REF-READS      a             b
//            //    MUT-READS      c             d
//            int a = 0, b = 0, c = 0, d = 0;
//            for (byte base : otherAlleles) {
//                if (base == altAllele1) {
//                    a++;
//                } else if (base == altAllele2) {
//                    b++;
//                }
//            }
//            for (byte base : mutantAlleles) {
//                if (base == altAllele1) {
//                    c++;
//                } else if (base == altAllele2) {
//                    d++;
//                }
//            }
//
//            // calculate the maximum p-value we could get
//            //               OFFSET-ALT1   OFFSET-ALT2
//            //    REF-READS      a+b           0
//            //    MUT-READS       0           c+d
//            final double mostSignificantP = fisher.getTwoTailedP(a+b,0,0,c+d);
//
//            if (mostSignificantP < SKEW_FISHER_PVALUE_CUTOFF) {
//                final double p = fisher.getTwoTailedP(a,b,c,d);
//                skewLodOffsets.put(offset, p);
////                System.out.println( "Offset: " + offset + " msp: " + mostSignificantP + " (a,b,c,d): (" + a + "," + b + "," + c + "," + d + ") p: "+ p);
//            }
//
//        }
//
//        // perform Bonferroni Correction
//        int n = skewLodOffsets.size();
//        //System.out.println("Bonferroni n:"+n);
//        for(Integer offset : skewLodOffsets.keySet()) {
//            double p = skewLodOffsets.get(offset);
//            skewLodOffsets.put(offset, p * n);
//        }
//
//        return skewLodOffsets;
//
//    }

    int MAX_READ_MISMATCH_QUALITY_SCORE_SUM = 100;
    // TODO: make this parameterizable
    private static Character MAPPED_BY_MATE = 'M';
    private LocusReadPile filterReads(final ReferenceContext ref, final LocusReadPile pile, final char altAllele, final String reference, final long leftmostIndex) {
        ArrayList<PileupElement> newPileupElements = new ArrayList<PileupElement>();

        for ( PileupElement p : pile.finalPileup ) {
            final GATKSAMRecord read = p.getRead();
            final int offset = p.getOffset();

            int mismatchQualitySum =
                    AlignmentUtils.mismatchesInRefWindow(p, ref, false, true);

            // do we have to many mismatches overall?
            if (mismatchQualitySum > this.MAX_READ_MISMATCH_QUALITY_SCORE_SUM) {
                //if (read.getReadString().charAt(offset) == altAllele) { out.println("MAXERR " + read.getReadName() + " that supported altAllele with " + mismatchScore + " mmqsum"); }
                continue;
            }

            // is this a heavily clipped read?
            if (isReadHeavilyClipped(read)) {
                continue;
            }

            // was this read ONLY placed because it's mate was uniquely placed? (supplied by BWA)
            if (MAPPED_BY_MATE.equals(read.getAttribute("XT"))) {
                continue;
            }
            // if we're here... we passed all the read filters!
            newPileupElements.add(new PileupElement(read, p.getOffset()));


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

    private int mismatchQualitySum(final SAMRecord read, final String refSeq, final int refIndex) {
        return mismatchQualitySum(read, refSeq, refIndex, 0);
    }

    private int mismatchQualitySum(final SAMRecord read, final String refSeq, final int refIndex, final int minMismatchQualityScore) {
        final List<Mismatch> mismatches = getMismatches(read, refSeq, refIndex);
        return mismatchQualitySum(mismatches, minMismatchQualityScore);
    }

    private int mismatchQualitySum(final List<Mismatch> mismatches) {
        return mismatchQualitySum(mismatches, 0);
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
                        // FIXME: what is this case????
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
