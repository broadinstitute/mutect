package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import net.sf.picard.util.FormatUtil;
import net.sf.samtools.util.StringUtil;

import java.util.ArrayList;
import java.util.List;

/**
 *  Encapsulates the generation of a CallStats line based on the information in the Candidate
 */
public class CallStatsGenerator {
    private static final String TAB = "\t";
    private final FormatUtil fmt = new FormatUtil();

    private int[] minimalHeaderIndicies;

    private static final String[] COMPLETE_CALL_STATS_HEADER =
            new String[]{
                    "contig","position","context","ref_allele","alt_allele","tumor_name","normal_name","score","dbsnp_site",
                    "covered", "power", "tumor_power", "normal_power", "normal_power_nsp", "normal_power_wsp",
                    "total_pairs","improper_pairs","map_Q0_reads",
                    "init_t_lod","t_lod_fstar", "t_lod_fstar_forward", "t_lod_fstar_reverse", "tumor_f", "contaminant_fraction","contaminant_lod","t_q20_count", "t_ref_count","t_alt_count","t_ref_sum","t_alt_sum","t_ref_max_mapq","t_alt_max_mapq","t_ins_count","t_del_count",
                    "normal_best_gt","init_n_lod","normal_f",
                    "n_q20_count", "n_ref_count","n_alt_count","n_ref_sum","n_alt_sum",
                    "power_to_detect_positive_strand_artifact", "power_to_detect_negative_strand_artifact",
                    "strand_bias_counts",
                    "tumor_alt_fpir_median", "tumor_alt_fpir_mad","tumor_alt_rpir_median","tumor_alt_rpir_mad",
                    "observed_in_normals_count", "failure_reasons","judgement"
            };

    private static final String[] MINIMAL_CALL_STATS_HEADER =
            new String[]{
                    "contig","position","context","ref_allele","alt_allele","tumor_name","normal_name","score","dbsnp_site",
                    "covered", "power", "tumor_power", "normal_power",
                    "total_pairs","improper_pairs","map_Q0_reads",
                    "t_lod_fstar","tumor_f","contaminant_fraction","contaminant_lod",
                    "t_ref_count","t_alt_count","t_ref_sum","t_alt_sum","t_ref_max_mapq","t_alt_max_mapq","t_ins_count","t_del_count",
                    "normal_best_gt","init_n_lod", "n_ref_count","n_alt_count","n_ref_sum","n_alt_sum",
                    "observed_in_normals_count", "failure_reasons","judgement"
            };

    private boolean enableExtendedOutput;

    public CallStatsGenerator(boolean enableExtendedOutput) {
        this.enableExtendedOutput = enableExtendedOutput;
    }

    public String generateHeader() {
        String header;
        if (this.enableExtendedOutput) {
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
        return header;
    }

    public String generateCallStats(CandidateMutation candidate) {
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
        if (candidate.isDbsnpSite() && candidate.isCosmicSite()) {
            siteInfo = "DBSNP+COSMIC";
        }

        StringBuilder sb = new StringBuilder();
        int[] ci = candidate.getStrandContingencyTable();
        sb.append("(");
        sb.append(ci[0]).append(",");
        sb.append(ci[1]).append(",");
        sb.append(ci[2]).append(",");
        sb.append(ci[3]).append(")");
        String strandInfo = sb.toString();

        String[] msg = new String[] {
                candidate.getLocation().getContig(),
                format(candidate.getLocation().getStart()),
                candidate.getSequenceContext(),
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
                format(candidate.getNormalPowerNoSNPPrior()),
                format(candidate.getNormalPowerWithSNPPrior()),
                format(candidate.getTotalPairs()),
                format(candidate.getImproperPairs()),
                format(candidate.getMapQ0Reads()),
                format(candidate.getInitialTumorLod()),
                format(candidate.getTumorLodFStar()),
                format(candidate.getTumorLodFStarForward()),
                format(candidate.getTumorLodFStarReverse()),
                format(candidate.getTumorF()),
                format(candidate.getContaminationFraction()),
                format(candidate.getContaminantLod()),
                format(candidate.getTumorQ20Count()),
                format(candidate.getInitialTumorRefCounts()),
                format(candidate.getInitialTumorAltCounts()),
                format(candidate.getInitialTumorRefQualitySum()),
                format(candidate.getInitialTumorAltQualitySum()),
                format(candidate.getTumorRefMaxMapQ()),
                format(candidate.getTumorAltMaxMapQ()),
                format(candidate.getTumorInsertionCount()),
                format(candidate.getTumorDeletionCount()),
                format(candidate.getInitialNormalBestGenotype().toString()),
                format(candidate.getInitialNormalLod()),
                format(candidate.getNormalF()),
                format(candidate.getNormalQ20Count()),
                format(candidate.getInitialNormalRefCounts()),
                format(candidate.getInitialNormalAltCounts()),
                format(candidate.getInitialNormalRefQualitySum()),
                format(candidate.getInitialNormalAltQualitySum()),
                format(candidate.getPowerToDetectPositiveStrandArtifact()),
                format(candidate.getPowerToDetectNegativeStrandArtifact()),
                format(strandInfo),
                candidate.getTumorForwardOffsetsInReadMedian()==null?"n/a":format(candidate.getTumorForwardOffsetsInReadMedian()),
                candidate.getTumorForwardOffsetsInReadMad()==null?"n/a":format(candidate.getTumorForwardOffsetsInReadMad()),
                candidate.getTumorReverseOffsetsInReadMedian()==null?"n/a":format(candidate.getTumorReverseOffsetsInReadMedian()),
                candidate.getTumorReverseOffsetsInReadMad()==null?"n/a":format(candidate.getTumorReverseOffsetsInReadMad()),
                format(candidate.getCountOfNormalsObservedIn()),
                StringUtil.join(",", candidate.getRejectionReasons().toArray(new String[candidate.getRejectionReasons().size()])),
                keepString
        };

        if (this.enableExtendedOutput) {
            for(int i=0; i<msg.length; i++) {
                if (msg[i]==null) { System.out.println("ERROR: found null in output at position " + i); }
            }
            return StringUtil.join(TAB, msg);
        } else {
            List<String> output = new ArrayList<String>();
            for(int index : minimalHeaderIndicies) {
                output.add(msg[index]);
            }
            return StringUtil.join(TAB, output);
        }
    }

    private String format(String s) { return s; }
    private String format(Integer i) { return fmt.format(i); }
    private String format(Double d) {
        if (d == null) { return "n/a"; }

        String s = fmt.format(d);
        return ("-0".equals(s))?"0":s;
    }
}
