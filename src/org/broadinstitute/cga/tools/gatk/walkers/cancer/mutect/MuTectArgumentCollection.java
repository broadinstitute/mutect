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
