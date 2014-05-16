/*
 * By downloading the PROGRAM you agree to the following terms of use:
 *
 * BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
 * FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
 *
 * This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
 * WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
 * WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
 * NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
 *
 * 1. DEFINITIONS
 * 1.1	"PROGRAM" shall mean copyright in the object code and source code known as MuTect and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.or/cancer/cga/mutect on the EFFECTIVE DATE.
 *
 * 2. LICENSE
 * 2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.
 * LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
 * The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
 * 2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
 * 2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
 *
 * 3. OWNERSHIP OF INTELLECTUAL PROPERTY
 * LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
 *
 * Copyright 2012 Broad Institute, Inc.
 * Notice of attribution:  The MuTect program was made available through the generosity of the Cancer Genome Analysis group at the Broad Institute, Inc. and is published at doi: 10.1038/nbt.2514.
 *
 * LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
 *
 * 4. INDEMNIFICATION
 * LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, ("Indemnitees"), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
 *
 * 5. NO REPRESENTATIONS OR WARRANTIES
 * THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
 * IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 *
 * 6. ASSIGNMENT
 * This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
 *
 * 7. MISCELLANEOUS
 * 7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
 * 7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
 * 7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
 * 7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
 * 7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
 * 7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
 * 7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
 */

package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMRecord;
import org.apache.commons.lang.math.NumberUtils;
import org.broadinstitute.cga.tools.gatk.utils.CGAAlignmentUtils;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.text.TextFormattingUtils;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;

import java.io.PrintStream;
import java.util.*;

@PartitionBy(PartitionType.LOCUS)
@Reference(window=@Window(start=-1* MuTect.REFERENCE_HALF_WINDOW_LENGTH,stop= MuTect.REFERENCE_HALF_WINDOW_LENGTH))
@By(DataSource.REFERENCE)
public class MuTect extends LocusWalker<Integer, Integer>  {
    public enum SampleType {TUMOR, NORMAL}
    public static final int REFERENCE_HALF_WINDOW_LENGTH = 150;
    public static final String BAM_TAG_TUMOR = "tumor";
    public static final String BAM_TAG_NORMAL = "normal";


    @ArgumentCollection private MuTectArgumentCollection MTAC = new MuTectArgumentCollection();

    /***************************************/
    // Call-stats output
    /***************************************/
    @Output(doc="Call-stats output")
    PrintStream out;

    @Output(doc="VCF output of mutation candidates",shortName="vcf", fullName="vcf", required=false, defaultToStdout = false)
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
    @Output(fullName="coverage_file", shortName="cov", doc="write out coverage in WIGGLE format to this file", required=false, defaultToStdout=false)
    public PrintStream COVERAGE_FILE = null;

    @Output(fullName="coverage_20_q20_file", shortName="cov_q20", doc="write out 20x of Q20 coverage in WIGGLE format to this file", required=false, defaultToStdout=false)
    public PrintStream COVERAGE_20_Q20_FILE = null;

    @Output(fullName="power_file", shortName="pow", doc="write out power in WIGGLE format to this file", required=false, defaultToStdout=false)
    public PrintStream POWER_FILE = null;

    @Output(fullName="tumor_depth_file", shortName="tdf", doc="write out tumor read depth in WIGGLE format to this file", required=false, defaultToStdout=false)
    public PrintStream TUMOR_DEPTH_FILE = null;

    @Output(fullName="normal_depth_file", shortName="ndf", doc="write out normal read depth in WIGGLE format to this file", required=false, defaultToStdout=false)
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

    public static class PileupComparatorByAltRefQual implements Comparator<PileupElement> {
        private byte alt;

        public PileupComparatorByAltRefQual(byte alt) {
            this.alt = alt;
        }

        public int compare(PileupElement o1, PileupElement o2) {
            return internalCompare(o1.getBase(), o1.getQual(), o2.getBase(), o2.getQual());
        }

        public int internalCompare(byte base1, byte qual1, byte base2, byte qual2) {
            // if the bases are the same, the higher quality score comes first
            if (base1 == base2) {
                if (qual1 == qual2) { return 0; }
                return (qual1 > qual2)?-1:1;

            // if the bases are not the same, then the alternate is first
            } else {
                if (base1 == alt) {
                    return -1;
                } else if (base2 == alt) {
                    return 1;
                } else {
                    return base1 < base2?-1:1;
                }

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
    public boolean includeReadsWithDeletionAtLoci() {
        return true;
    }

    @Override
    public void initialize() {
        if (MTAC.NOOP) {
            return;
        }

        //setting version info
        // TODO: refactor into getMuTectVersion()
        final String gatkVersion = CommandLineGATK.getVersionNumber();
        ResourceBundle resources = TextFormattingUtils.loadResourceBundle("MuTectText");
        final String mutectVersion = resources.containsKey("version")? resources.getString("version") : "<unknown>";        
        final String combinedVersion = "MuTect:"+mutectVersion+" Gatk:"+gatkVersion;

        logger.info("VERSION INFO: " + combinedVersion);

        refReader = this.getToolkit().getReferenceDataSource().getReference();
        callStatsGenerator = new CallStatsGenerator(MTAC.ENABLE_QSCORE_OUTPUT);

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
        out.println("##"+combinedVersion);
        out.println(callStatsGenerator.generateHeader());

        // initialize the VCF output
        if (vcf != null) {
            // TODO: fix for multisample mode
            Set<String> samples = new HashSet<String>();
            samples.add(MTAC.TUMOR_SAMPLE_NAME);
            samples.add(MTAC.NORMAL_SAMPLE_NAME);
            Set<VCFHeaderLine> headerInfo = VCFGenerator.getVCFHeaderInfo();
            vcf.writeHeader(new VCFHeader(headerInfo, samples));
        }

        lastTime = System.currentTimeMillis();
    }

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

        // get sequence context around mutation
        String sequenceContext = SequenceUtils.createSequenceContext(ref, 3);

        // only process bases where the reference is [ACGT], because the FASTA for HG18 has N,M and R!
        final char upRef = Character.toUpperCase(ref.getBaseAsChar());
        if (upRef != 'A' && upRef != 'C' && upRef != 'G' && upRef != 'T') {
            return -1;
        }

        try {

            Map<SampleType, ReadBackedPileup> pileupMap = getPileupsBySampleType(pileup);



            final LocusReadPile tumorReadPile =  new LocusReadPile(pileupMap.get(SampleType.TUMOR), upRef, MTAC.MIN_QSCORE, MIN_QSUM_QSCORE, false, MTAC.ARTIFACT_DETECTION_MODE, MTAC.ENABLE_QSCORE_OUTPUT);
            final LocusReadPile normalReadPile = new LocusReadPile(pileupMap.get(SampleType.NORMAL), upRef, MTAC.MIN_QSCORE, 0, this.USE_MAPQ0_IN_NORMAL_QSCORE, true, MTAC.ENABLE_QSCORE_OUTPUT);

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

            int totalReads =
                    tumorReadPile.qualityScoreFilteredPileup.depthOfCoverage() +
                            normalReadPile.qualityScoreFilteredPileup.depthOfCoverage();

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
                candidate.setMapQ0Reads(mapQ0Reads);
                candidate.setTotalReads(totalReads);
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

                Collections.sort(peList, new PileupComparatorByAltRefQual((byte)altAllele));
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

                candidate.setNormalAltQualityScores(normQs.getBaseQualityScores(altAllele));
                candidate.setNormalRefQualityScores(normQs.getBaseQualityScores(upRef));

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

                candidate.setTumorAltQualityScores(t2.qualitySums.getBaseQualityScores(altAllele));
                candidate.setTumorRefQualityScores(t2.qualitySums.getBaseQualityScores(upRef));

                VariableAllelicRatioGenotypeLikelihoods t2Gl = t2.calculateLikelihoods(t2.finalPileup);
                candidate.setInitialTumorLod(t2.getAltVsRef(t2Gl, upRef, altAllele));
                candidate.setInitialTumorReadDepth(t2.finalPileupReads.size());

                candidate.setTumorF(t2.estimateAlleleFraction(upRef, altAllele));
                double tumorLod2 = t2.calculateAltVsRefLOD((byte)altAllele, candidate.getTumorF(), 0);
                candidate.setTumorLodFStar(tumorLod2);

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

                candidate.setStrandContingencyTable(SequenceUtils.getStrandContingencyTable(forwardPileup, reversePileup, (byte) upRef, (byte) altAllele));

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
                final LocusReadPile mutantPile = new LocusReadPile(mutantPileup, altAllele, 0, 0, MTAC.ENABLE_QSCORE_OUTPUT);
                final LocusReadPile refPile =  new LocusReadPile(referencePileup, altAllele, 0, 0, MTAC.ENABLE_QSCORE_OUTPUT);

                // Set the maximum observed mapping quality score for the reference and alternate alleles
                int[] rmq = referencePileup.getMappingQuals();
                candidate.setTumorRefMaxMapQ((rmq.length==0)?0:NumberUtils.max(rmq));

                int[] amq = mutantPileup.getMappingQuals();
                candidate.setTumorAltMaxMapQ((amq.length==0)?0:NumberUtils.max(amq));


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

        if (candidate.getTotalReads() > 0 && ((float)candidate.getMapQ0Reads() / (float)candidate.getTotalReads()) >= MTAC.FRACTION_MAPQ0_THRESHOLD) {
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


    protected Map<SampleType, ReadBackedPileup> getPileupsBySampleType(ReadBackedPileup pileup) {
        Map<SampleType, ReadBackedPileup> result = new HashMap<SampleType, ReadBackedPileup>();

        ArrayList<PileupElement> tumorPileupElements = new ArrayList<PileupElement>();
        ArrayList<PileupElement> normalPileupElements = new ArrayList<PileupElement>();

        for (PileupElement p : pileup ) {
            final byte base = p.getBase();

            if (base == ((byte)'N') || base == ((byte)'n')) {
                continue;
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
                new ReadBackedPileupImpl(pileup.getLocation(), normalPileupElements);

        ReadBackedPileup tumorPileup =
                new ReadBackedPileupImpl(pileup.getLocation(), tumorPileupElements);

        // if the user specified a single tumor sample name, only look at that sample
        if (MTAC.BAM_TUMOR_SAMPLE_NAME != null) {
            tumorPileup = tumorPileup.getPileupForSample(MTAC.BAM_TUMOR_SAMPLE_NAME);
        }

        result.put(SampleType.NORMAL, normalPileup);
        result.put(SampleType.TUMOR, tumorPileup);

        return result;
    }

    int MAX_READ_MISMATCH_QUALITY_SCORE_SUM = 100;
    private static Character MAPPED_BY_MATE = 'M';
    IndexedFastaSequenceFile refReader;

    private LocusReadPile filterReads(final ReferenceContext ref, final ReadBackedPileup pile, boolean filterMateRescueReads) {
        ArrayList<PileupElement> newPileupElements = new ArrayList<PileupElement>();

        for ( PileupElement p : pile ) {
            final GATKSAMRecord read = p.getRead();

            int mismatchQualitySum =
                    CGAAlignmentUtils.mismatchesInRefWindow(p, ref, false, true);

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
			            newPileupElements.add(new PileupElement(p));


        }
        ReadBackedPileup newPileup =
                new ReadBackedPileupImpl(ref.getLocus(), newPileupElements);


        return new LocusReadPile(newPileup, (char)ref.getBase(), 0, 0, MTAC.ENABLE_QSCORE_OUTPUT);
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
