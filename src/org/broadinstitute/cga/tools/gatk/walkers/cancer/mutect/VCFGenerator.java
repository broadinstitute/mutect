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

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class VCFGenerator {
    public static Set<VCFHeaderLine> getVCFHeaderInfo() {
        Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();

        headerInfo.add(new VCFFilterHeaderLine("triallelic_site", "Triallelic sites generate false positives, and are currently filtered out by default"));
        headerInfo.add(new VCFFilterHeaderLine("fstar_tumor_lod", ""));
        headerInfo.add(new VCFFilterHeaderLine("nearby_gap_events", "Nearby misaligned small insertion and deletion events"));
        headerInfo.add(new VCFFilterHeaderLine("possible_contamination", ""));
        headerInfo.add(new VCFFilterHeaderLine("germline_risk", ""));
        headerInfo.add(new VCFFilterHeaderLine("normal_lod", ""));
        headerInfo.add(new VCFFilterHeaderLine("alt_allele_in_normal", "ALT allele observed in control"));
        headerInfo.add(new VCFFilterHeaderLine("clustered_read_position", ""));
        headerInfo.add(new VCFFilterHeaderLine("strand_artifact", "Majority of ALT alleles are observed in a single direction of reads"));
        headerInfo.add(new VCFFilterHeaderLine("poor_mapping_region_mapq0", "≥50% reads in tumor and normal have a mapping quality of zero"));
        headerInfo.add(new VCFFilterHeaderLine("poor_mapping_region_alternate_allele_mapq", "No ALT allele reads with mapping quality score ≥ 20"));
        headerInfo.add(new VCFFilterHeaderLine("seen_in_panel_of_normals", "Present in two or more normal samples"));
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
                        .attribute(VCFConstants.RMS_BASE_QUALITY_KEY, "." );    // TODO: confirm TCGA-compliancy

        VariantContextBuilder vc =
                new VariantContextBuilder("", l.getContig(), l.getStart(), l.getStop(), alleles);

        vc.filters( m.isRejected() ? new HashSet<String>(m.getRejectionReasons()) : new HashSet<String>(Arrays.asList("PASS")));
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
