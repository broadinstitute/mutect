package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.cancer.MuTect
import java.io.File
import org.broadinstitute.sting.queue.extensions.gatk.ArgumentDefinitionField.InputTaggedFileDefinitionField
import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.commandline.ClassType
import org.broadinstitute.sting.commandline

class MuTectPipeline extends QScript {

  @Input(doc="input tumor BAM file - or list of BAM files", shortName="tb", required=true)
  var tumor_bams: Seq[File] = Seq()

  @Input(doc="input normal BAM file - or list of BAM files", shortName="nb", required=true)
  var normal_bams: Seq[File] = Seq()

  @Argument(doc="Output head", shortName="o", required=true)
  var output_head: String = _ 

  // ---------------------------------- optional arguments ----------------------------------
  @Input(doc="Reference fasta file", fullName="reference", shortName="R", required=false)
  var reference: File = _

  @Input(doc="DBSNP or known callset to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=false)
  var dbSNP: Seq[File] = Seq()

  @Input(doc="Panel Of Normals or known artifact sites to use (must be in VCF format)", fullName="panel_of_normals", shortName="pon", required=false)
  var pon: Seq[File] = Seq()

  @Input(doc="COSMIC sites to use (must be in VCF format)", fullName="cosmic", shortName="C", required=false)
  var cosmic: Seq[File] = Seq()

  @Argument(doc="Configuration (hg18-wex|hg19-wex)", shortName="config", required=false)
  var config: String = ""

  @Argument(doc="force caller output  ", shortName="fo", required = false)
  var forceOutput: Boolean = false

  @Input(doc="intervals to use for calling", required = false)
  var intervals: Seq[File] = Nil

  @Argument(doc="scatter count", required = false)
  var sg = 1

  @Argument(doc="tumor sample name", required = false)
  var tumor_name = output_head + "-Tumor"
  
  @Argument(doc="normal sample name", required = false)
  var normal_name = output_head + "-Normal"

  @ClassType(classOf[Float])
  @Argument(doc="Contamination Fraction", shortName="c", required = false)
  var contamination = Some(0.02f)



  def script = {

  if (config == "hg19-wex") {
    reference = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta")
    pon :+= new File("/xchip/cga/reference/hg19/refseq_exome_10bp_hg19_300_1kg_normal_panel.vcf")
    dbSNP :+= new File("/xchip/cga/reference/hg19/dbsnp_134_b37.leftAligned.vcf")
    cosmic :+= new File("/xchip/cga/reference/hg19/hg19_cosmic_v54_120711.vcf")
    intervals :+= new File("/xchip/cga/reference/hg19/gaf_20111020+broad_wex_1.1_hg19.bed")
  } else if (config == "hg18-wex") {
    reference = new File("/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta")
    pon :+= new File("/xchip/cga/reference/hg18/refseq_exome_10bp_hg18_300_1kg_normal_panel.vcf")
    dbSNP :+= new File("/xchip/cga/reference/hg18/dbsnp_132.hg18.vcf")
    cosmic :+= new File("/xchip/cga/reference/hg18/hg18_cosmic_v54_120711.vcf")
    intervals :+= new File("/xchip/cga1/kcibul/analysis/exome_targets/whole_exome_agilent_1.1_refseq_plus_3_boosters_plus_tcga_6k_plus_10bp_padding_minus_mito.Homo_sapiens_assembly18.targets.interval_list")
  }

	val m = new MuTect()
  m.normal_panel = pon
  m.cosmic = cosmic
  m.dbsnp = dbSNP
	m.enable_extended_output = true
  m.force_output = forceOutput
	m.intervals = intervals
	m.reference_sequence = reference
  m.fraction_contamination = contamination
	m.tumor_sample_name = tumor_name
	m.normal_sample_name = normal_name
	m.out = output_head + ".call_stats.txt"
  m.coverage_file = output_head + ".coverage.wig.txt"
  m.tumor_depth_file = output_head + ".tumor.depth.wig.txt"
  m.normal_depth_file = output_head + ".normal.depth.wig.txt"

	for (bam <- tumor_bams) {
		m.input_file :+= new TaggedFile(bam, "tumor")
	}
	for (bam <- normal_bams) {
		m.input_file :+= new TaggedFile(bam, "normal")
	}

  m.scatterCount = sg
	m.memoryLimit = 4

	
	add(m)

  }

}
