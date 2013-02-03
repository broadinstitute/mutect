package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import org.broadinstitute.sting.commandline.*;

public class MuTectArgumentCollection {
    @Hidden
    @Argument(fullName = "noop", required = false, doc="used for debugging, basically exit as soon as we get the reads")
    public boolean NOOP = false;

    @Argument(fullName = "enable_extended_output", required = false, doc="add many additional columns of statistics to the output file")
    public boolean ENABLE_EXTENDED_OUTPUT = false;

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

    @Argument(fullName = "only_passing_calls", required = false, doc="only emit passing calls")
    public boolean ONLY_PASSING_CALLS = false;

    @Argument(fullName = "initial_tumor_lod", required = false, doc = "Initial LOD threshold for calling tumor variant")
    public float INITIAL_TUMOR_LOD_THRESHOLD = 4.0f;

    @Argument(fullName = "tumor_lod", required = false, doc = "LOD threshold for calling tumor variant")
    public float TUMOR_LOD_THRESHOLD = 6.3f;

    @Argument(fullName = "fraction_contamination", required = false, doc = "estimate of fraction (0-1) of physical contamination with other unrelated samples")
    public float FRACTION_CONTAMINATION = 0.02f;

    @Argument(fullName = "minimum_mutation_cell_fraction", required = false,
            doc = "minimum fraction of cells which are presumed to have a mutation, used to handle non-clonality and contamination")
    public float MINIMUM_MUTATION_CELL_FRACTION = 0.00f;

    @Argument(fullName = "normal_lod", required = false, doc = "LOD threshold for calling normal non-germline")
    public float NORMAL_LOD_THRESHOLD = 2.2f;

    @Hidden
    @Argument(fullName = "normal_artifact_lod", required = false, doc = "LOD threshold for calling normal non-variant")
    public float NORMAL_ARTIFACT_LOD_THRESHOLD = 1.0f;

    @Hidden
    @Argument(fullName = "strand_artifact_lod", required = false, doc = "LOD threshold for calling strand bias")
    public float STRAND_ARTIFACT_LOD_THRESHOLD = 2.0f;

    @Hidden
    @Argument(fullName = "strand_artifact_power_threshold", required = false, doc = "power threshold for calling strand bias")
    public float STRAND_ARTIFACT_POWER_THRESHOLD = 0.9f;

    @Argument(fullName = "dbsnp_normal_lod", required = false, doc = "LOD threshold for calling normal non-variant at dbsnp sites")
    public float NORMAL_DBSNP_LOD_THRESHOLD = 5.5f;

    @Argument(fullName = "minimum_normal_allele_fraction", required = false, doc = "minimum allele fraction to be considered in normal, useful for normal sample contaminated with tumor")
    public float MINIMUM_NORMAL_ALLELE_FRACTION = 0.00f;

    @Argument(fullName = "tumor_f_pretest", required = false, doc = "for computational efficiency, reject sites with allelic fraction below this threshold")
    public float TUMOR_F_PRETEST = 0.005f;

    @Argument(fullName = "min_qscore", required = false, doc = "threshold for minimum base quality score")
    public int MIN_QSCORE = 5;

    @Argument(fullName = "gap_events_threshold", required = false, doc = "how many gapped events (ins/del) are allowed in proximity to this candidate")
    public int GAP_EVENTS_THRESHOLD = 3;

    @Argument(fullName = "heavily_clipped_read_fraction", required = false, doc = "if this fraction or more of the bases in a read are soft/hard clipped, do not use this read for mutation calling")
    public float HEAVILY_CLIPPED_READ_FRACTION = 0.30f;

    @Argument(fullName = "fraction_mapq0_threshold", required = false, doc = "threshold for determining if there is relatedness between the alt and ref allele read piles")
    public float FRACTION_MAPQ0_THRESHOLD = 0.5f;

    @Argument(fullName = "pir_median_threshold", required = false, doc="threshold for clustered read position artifact median")
    public double PIR_MEDIAN_THRESHOLD = 10;

    @Argument(fullName = "pir_mad_threshold", required = false, doc="threshold for clustered read position artifact MAD")
    public double PIR_MAD_THRESHOLD = 3;

    @Argument(fullName = "required_maximum_alt_allele_mapping_quality_score", required = false, doc="required minimum value for tumor alt allele maximum mapping quality score")
    public int REQUIRED_MAXIMUM_ALT_ALLELE_MAPPING_QUALITY_SCORE = 20;

    /** Parameters for ALT ALLELE IN NORMAL filter **/
    @Argument(fullName = "max_alt_alleles_in_normal_count", required = false, doc="threshold for maximum alternate allele counts in normal")
    public int MAX_ALT_ALLELES_IN_NORMAL_COUNT = 2;

    @Argument(fullName = "max_alt_alleles_in_normal_qscore_sum", required = false, doc="threshold for maximum alternate allele quality score sum in normal")
    public int MAX_ALT_ALLELES_IN_NORMAL_QSCORE_SUM = 20;

    @Argument(fullName = "max_alt_allele_in_normal_fraction", required = false, doc="threshold for maximum alternate allele fraction in normal")
    public double MAX_ALT_ALLELE_IN_NORMAL_FRACTION = 0.03;

    @Argument(fullName="power_constant_qscore", doc="Phred scale quality score constant to use in power calculations", required=false)
    public int POWER_CONSTANT_QSCORE = 30;

    @Argument(fullName="power_constant_af", doc="Allelic fraction constant to use in power calculations", required=false)
    public double POWER_CONSTANT_AF = 0.3f;

}
