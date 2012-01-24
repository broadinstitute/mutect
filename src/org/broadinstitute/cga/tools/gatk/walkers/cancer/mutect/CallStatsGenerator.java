package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import net.sf.picard.util.FormatUtil;
import net.sf.samtools.util.StringUtil;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 *  Encapsulates the generation of a CallStats line based on the information in the Candidate
 */
public class CallStatsGenerator {
    private static final String TAB = "\t";
    private final FormatUtil fmt = new FormatUtil();

    private int[] minimalHeaderIndicies;

    private static final String[] COMPLETE_CALL_STATS_HEADER =
            new String[]{
                "contig","position","ref_allele","alt_allele","tumor_name","normal_name","score","dbsnp_site",
                "covered", "power", "tumor_power", "normal_power",
                "total_pairs","improper_pairs","map_Q0_reads",
                "init_t_lod","t_lod_fstar","t_lod_fstar_forward", "t_lod_fstar_reverse", "tumor_f","contaminant_fraction","contaminant_lod","minimum_tumor_f", "t_q20_count", "t_ref_count","t_alt_count","t_ref_sum","t_alt_sum","t_ins_count","t_del_count",
                "normal_best_gt","init_n_lod","n_lod_fstar","normal_f","normal_artifact_lod","n_q20_count", "n_ref_count","n_alt_count","n_ref_sum","n_alt_sum",
                "at_risk_positive_direction_artifact", "at_risk_negative_direction_artifact", "powered_positive_direction_artifact", "power_to_detect_positive_strand_artifact", "powered_negative_direction_artifact", "power_to_detect_negative_strand_artifact",
                "perfect_strand_bias","strand_bias_counts","strand_bias",
                "classic_max_skew_lod", "classic_max_skew_lod_offset", "fisher_min_skew_pvalue", "fisher_min_skew_pvalue_offset",
                "tumor_qsrst_ms", "tumor_qsrst_pval", "tumor_rprst_ms", "tumor_rprst_pval",
                "tumor_alt_fpir_median", "tumor_alt_fpir_mad","tumor_alt_rpir_median","tumor_alt_rpir_mad","alt_fpir","alt_rpir",
                "powered_filters", "normal_artifact_power",
                "observed_in_normals_count", "failure_reasons","judgement"
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
        if (candidate.isDbsnpSite() && candidate.isCosmicSite()) {
            siteInfo = "DBSNP+COSMIC";
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
                        format(candidate.getContaminationFraction()),
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
                        format(candidate.getNormalArtifactLod()),
                        format(candidate.getNormalQ20Count()),
                        format(candidate.getInitialNormalRefCounts()),
                        format(candidate.getInitialNormalAltCounts()),
                        format(candidate.getInitialNormalRefQualitySum()),
                        format(candidate.getInitialNormalAltQualitySum()),
                        format(candidate.isPositiveDirectionAtRisk()?1:0),
                        format(candidate.isNegativeDirectionAtRisk()?1:0),
                        format(candidate.isPositiveDirectionPowered()?1:0),
                        format(candidate.getPowerToDetectPositiveStrandArtifact()),
                        format(candidate.isNegativeDirectionPowered()?1:0),
                        format(candidate.getPowerToDetectNegativeStrandArtifact()),
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
                        (candidate.getTumorAltForwardOffsetsInRead().size() > 500 ? "too_many": StringUtil.join(",", listToString(candidate.getTumorAltForwardOffsetsInRead()))),
                        (candidate.getTumorAltReverseOffsetsInRead().size() > 500 ? "too_many":StringUtil.join(",", listToString(candidate.getTumorAltReverseOffsetsInRead()))),
                        StringUtil.join(",", candidate.getPoweredFilters().toArray(new String[]{})),
                        format(candidate.getNormalArtifactPower()),
                        format(candidate.getCountOfNormalsObservedIn()),
                        StringUtil.join(",", candidate.getRejectionReasons().toArray(new String[]{})),
                        keepString
            };

        if (this.enableExtendedOutput) {
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

    private String[] listToString(List<Integer> ints) {
        String[] out = new String[ints.size()];
        for(int i=0; i<ints.size(); i++) {
            out[i] = ints.get(i).toString();
        }
        return out;
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


}
