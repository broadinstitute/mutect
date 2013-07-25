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

import net.sf.picard.util.FormatUtil;
import net.sf.samtools.util.StringUtil;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 *  Encapsulates the generation of a CallStats line based on the information in the Candidate
 */
public class CallStatsGenerator {
    private static final String TAB = "\t";
    private final FormatUtil fmt = new FormatUtil();

    private List<String> headers;

    public CallStatsGenerator(boolean enableQscoreOutput) {
        // set up the output columns (ordered!)
        this.headers = setupOutputColumns(enableQscoreOutput);
    }

    private List<String> setupOutputColumns(boolean enableQscoreOutput) {
        List<String> l = new ArrayList<String>();

        l.add("contig");
        l.add("position");
        l.add("context");
        l.add("ref_allele");
        l.add("alt_allele");
        l.add("tumor_name");
        l.add("normal_name");
        l.add("score");
        l.add("dbsnp_site");
        l.add("covered");
        l.add("power");
        l.add("tumor_power");
        l.add("normal_power");
        l.add("normal_power_nsp");
        l.add("normal_power_wsp");
        l.add("total_reads");
        l.add("map_Q0_reads");
        l.add("init_t_lod");
        l.add("t_lod_fstar");
        l.add("t_lod_fstar_forward");
        l.add("t_lod_fstar_reverse");
        l.add("tumor_f");
        l.add("contaminant_fraction");
        l.add("contaminant_lod");
        l.add("t_q20_count");
        l.add("t_ref_count");
        l.add("t_alt_count");
        l.add("t_ref_sum");
        l.add("t_alt_sum");
        l.add("t_ref_max_mapq");
        l.add("t_alt_max_mapq");
        l.add("t_ins_count");
        l.add("t_del_count");
        l.add("normal_best_gt");
        l.add("init_n_lod");
        l.add("normal_f");
        l.add("n_q20_count");
        l.add("n_ref_count");
        l.add("n_alt_count");
        l.add("n_ref_sum");
        l.add("n_alt_sum");
        l.add("power_to_detect_positive_strand_artifact");
        l.add("power_to_detect_negative_strand_artifact");
        l.add("strand_bias_counts");
        l.add("tumor_alt_fpir_median");
        l.add("tumor_alt_fpir_mad");
        l.add("tumor_alt_rpir_median");
        l.add("tumor_alt_rpir_mad");
        l.add("observed_in_normals_count");

        if (enableQscoreOutput) {
            l.add("tumor_ref_qscores");
            l.add("tumor_alt_qscores");
            l.add("normal_ref_qscores");
            l.add("normal_alt_qscores");
        }

        l.add("failure_reasons");
        l.add("judgement");
        return l;
    }

    public String generateHeader() {
        return StringUtil.join(TAB, this.headers);
    }

    public String generateCallStats(CandidateMutation candidate) {
        HashMap<String, String> d = new HashMap<String, String>();

        String keepString = "REJECT";
        if (!candidate.isRejected()) {
            keepString = "KEEP";
        }

        String siteInfo =
                getSiteInfoString(candidate.isDbsnpSite(), candidate.isCosmicSite());

        String strandInfo =
                getStrandTableString(candidate.getStrandContingencyTable());

        d.put("contig", candidate.getLocation().getContig());
        d.put("position", format(candidate.getLocation().getStart()));
        d.put("context", candidate.getSequenceContext());
        d.put("ref_allele", ""+candidate.getRefAllele());
        d.put("alt_allele", ""+candidate.getAltAllele());
        d.put("tumor_name", candidate.getTumorSampleName());
        d.put("normal_name", candidate.getNormalSampleName());
        d.put("score", format(candidate.getScore()));
        d.put("dbsnp_site", siteInfo);
        d.put("covered", (candidate.isCovered()?"COVERED":"UNCOVERED"));
        d.put("power", format(candidate.getPower()));
        d.put("tumor_power", format(candidate.getTumorPower()));
        d.put("normal_power", format(candidate.getNormalPower()));
        d.put("normal_power_nsp", format(candidate.getNormalPowerNoSNPPrior()));
        d.put("normal_power_wsp", format(candidate.getNormalPowerWithSNPPrior()));
        d.put("total_reads", format(candidate.getTotalReads()));
        d.put("map_Q0_reads", format(candidate.getMapQ0Reads()));
        d.put("init_t_lod", format(candidate.getInitialTumorLod()));
        d.put("t_lod_fstar", format(candidate.getTumorLodFStar()));
        d.put("t_lod_fstar_forward", format(candidate.getTumorLodFStarForward()));
        d.put("t_lod_fstar_reverse", format(candidate.getTumorLodFStarReverse()));
        d.put("tumor_f", format(candidate.getTumorF()));
        d.put("contaminant_fraction", format(candidate.getContaminationFraction()));
        d.put("contaminant_lod", format(candidate.getContaminantLod()));
        d.put("t_q20_count", format(candidate.getTumorQ20Count()));
        d.put("t_ref_count", format(candidate.getInitialTumorRefCounts()));
        d.put("t_alt_count", format(candidate.getInitialTumorAltCounts()));
        d.put("t_ref_sum", format(candidate.getInitialTumorRefQualitySum()));
        d.put("t_alt_sum", format(candidate.getInitialTumorAltQualitySum()));
        d.put("t_ref_max_mapq", format(candidate.getTumorRefMaxMapQ()));
        d.put("t_alt_max_mapq", format(candidate.getTumorAltMaxMapQ()));
        d.put("t_ins_count", format(candidate.getTumorInsertionCount()));
        d.put("t_del_count", format(candidate.getTumorDeletionCount()));
        d.put("normal_best_gt", format(candidate.getInitialNormalBestGenotype().toString()));
        d.put("init_n_lod", format(candidate.getInitialNormalLod()));
        d.put("normal_f", format(candidate.getNormalF()));
        d.put("n_q20_count", format(candidate.getNormalQ20Count()));
        d.put("n_ref_count", format(candidate.getInitialNormalRefCounts()));
        d.put("n_alt_count", format(candidate.getInitialNormalAltCounts()));
        d.put("n_ref_sum", format(candidate.getInitialNormalRefQualitySum()));
        d.put("n_alt_sum", format(candidate.getInitialNormalAltQualitySum()));
        d.put("power_to_detect_positive_strand_artifact", format(candidate.getPowerToDetectPositiveStrandArtifact()));
        d.put("power_to_detect_negative_strand_artifact", format(candidate.getPowerToDetectNegativeStrandArtifact()));
        d.put("strand_bias_counts", format(strandInfo));
        d.put("tumor_alt_fpir_median", candidate.getTumorForwardOffsetsInReadMedian()==null?"n/a":format(candidate.getTumorForwardOffsetsInReadMedian()));
        d.put("tumor_alt_fpir_mad", candidate.getTumorForwardOffsetsInReadMad()==null?"n/a":format(candidate.getTumorForwardOffsetsInReadMad()));
        d.put("tumor_alt_rpir_median", candidate.getTumorReverseOffsetsInReadMedian()==null?"n/a":format(candidate.getTumorReverseOffsetsInReadMedian()));
        d.put("tumor_alt_rpir_mad", candidate.getTumorReverseOffsetsInReadMad()==null?"n/a":format(candidate.getTumorReverseOffsetsInReadMad()));
        d.put("observed_in_normals_count", format(candidate.getCountOfNormalsObservedIn()));
        d.put("tumor_ref_qscores", format(candidate.getTumorRefQualityScores()));
        d.put("tumor_alt_qscores", format(candidate.getTumorAltQualityScores()));
        d.put("normal_ref_qscores", format(candidate.getNormalRefQualityScores()));
        d.put("normal_alt_qscores", format(candidate.getNormalAltQualityScores()));
        d.put("failure_reasons", StringUtil.join(",", candidate.getRejectionReasons().toArray(new String[candidate.getRejectionReasons().size()])));
        d.put("judgement", keepString);

        return generate(d);
    }

    private String generate(HashMap<String, String> d) {
        String[] msg = new String[headers.size()];

        for(int i=0; i<msg.length; i++) {
            String value = d.get(headers.get(i));
            if (value == null) { value = ""; }
            msg[i] = value;
        }
        return StringUtil.join(TAB, msg);
    }

    private String format(String s) { return s; }
    private String format(Integer i) { return fmt.format(i); }
    private String format(Double d) {
        if (d == null) { return "n/a"; }

        String s = fmt.format(d);
        return ("-0".equals(s))?"0":s;
    }

    private String format(List<Integer> ints) {
        if (ints == null || ints.size() == 0) {
            return "n/a";
        }

        return StringUtil.join(",", ints);
    }

    private String getSiteInfoString(boolean isDbsnpSite, boolean isCosmicSite) {
        String siteInfo = "NOVEL";
        if (isDbsnpSite) {
            siteInfo = "DBSNP";
        }
        if (isCosmicSite) {
            siteInfo = "COSMIC";
        }
        if (isDbsnpSite && isCosmicSite) {
            siteInfo = "DBSNP+COSMIC";
        }
        return siteInfo;
    }

    private String getStrandTableString(int[] ci) {
        StringBuilder sb = new StringBuilder();
        sb.append("(");
        sb.append(ci[0]).append(",");
        sb.append(ci[1]).append(",");
        sb.append(ci[2]).append(",");
        sb.append(ci[3]).append(")");
        return sb.toString();
    }

}
