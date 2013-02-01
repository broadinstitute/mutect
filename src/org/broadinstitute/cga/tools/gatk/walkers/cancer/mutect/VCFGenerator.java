package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.GenotypeBuilder;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class VCFGenerator {
    public static Set<VCFHeaderLine> getVCFHeaderInfo() {
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

    public static VariantContext generateVC(CandidateMutation m) {
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
            normalGenotype.attribute("SS", 0); // TODO: extract these TCGA specific attributes to a class
        }

        // add the genotype objects
        vc.genotypes(tumorGenotype.make(), normalGenotype.make());

        return vc.make();
    }
}
