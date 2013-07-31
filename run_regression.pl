#!/usr/bin/perl -w
use strict;

my $size = shift() || "micro";

my $REF="/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
my $DBSNP="/xchip/cga/reference/hg19/dbsnp_132.b37.vcf";
my $COSMIC="/xchip/cga/reference/hg19/hg19_cosmic_v54_120711.vcf";
my $PON="/xchip/cga/reference/hg19/refseq_exome_10bp_hg19_300_1kg_normal_panel.vcf";
my $MUTECT_JAR="../../gatk-protected/dist/GenomeAnalysisTK.jar";
my $OUTPUT_DIR="/tmp/regression_trial";

my $cmd =
			"mkdir -p $OUTPUT_DIR && java -Xmx1g " .
			"-jar $MUTECT_JAR " .
			"--analysis_type MuTect " .
			"--reference_sequence $REF " .
			"--dbsnp $DBSNP " .
			"--cosmic $COSMIC " .
			"--normal_panel $PON " .
			"--tumor_sample_name Tumor --normal_sample_name Normal " .
			"--fraction_contamination 0.01 " .
			"--force_output " .
			"-dt NONE " .
			"-I:normal testdata/HCC1143_BL.cghub.ccle.$size.bam " .
			"-I:tumor testdata/HCC1143.cghub.ccle.$size.bam " .
			"-o $OUTPUT_DIR/ccle.$size.call_stats.txt " .
			"-vcf $OUTPUT_DIR/ccle.$size.vcf " .
			"--tumor_depth_file $OUTPUT_DIR/ccle.$size.tumor.depth.wig.txt " .
			"--normal_depth_file $OUTPUT_DIR/ccle.$size.normal.depth.wig.txt " .
			"-L testdata/target.exons.$size.bed";

system($cmd) == 0 or die();
