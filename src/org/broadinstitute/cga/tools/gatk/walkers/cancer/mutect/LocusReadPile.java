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

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.genotyper.DiploidGenotype;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static java.lang.Math.pow;

public class LocusReadPile {
    protected ReadBackedPileup pileup;
    List<GATKSAMRecord> finalPileupReads;
    List<Integer> finalPileupOffsets;

    protected char refBase;
    protected int minQualityScore;
    protected int minQSumQualityScore;
    protected QualitySums qualitySums;

    public ReadBackedPileup initialPileup;
    protected ReadBackedPileup qualityScoreFilteredPileup;
    public ReadBackedPileup finalPileup;


    public ReadBackedPileup finalPileupPositiveStrand;
    public ReadBackedPileup finalPileupNegativeStrand;

    public static final int GAP_EVENT_PROXIMITY = 5; // a 11bp window
    protected int deletionsCount = 0;
    protected int insertionsCount = 0;

    public LocusReadPile(ReadBackedPileup pileup, char refBase, int minQualityScore, int minQSumQualityScore, boolean trackBaseQualityScores) {
        this(pileup, refBase, minQualityScore, minQSumQualityScore, false, false, trackBaseQualityScores);
    }

    public LocusReadPile(ReadBackedPileup pileup, char refBase, int minQualityScore, int minQSumQualityScore, boolean allowMapq0ForQualSum, boolean retainOverlapMismatches, boolean trackBaseQualityScores) {
        this.qualitySums = new QualitySums(trackBaseQualityScores);
        this.pileup = pileup;
        this.refBase = refBase;
        this.minQSumQualityScore = minQSumQualityScore;
        this.minQualityScore = minQualityScore;

        ReadBackedPileup noOverlapPileup;
        if (!retainOverlapMismatches) {
            noOverlapPileup = getOverlappingFragmentFilteredPileup(pileup, (byte) refBase);
        } else {
            noOverlapPileup = getOverlappingFragmentFilteredPileupButPreferMismatches(pileup, (byte) refBase);
        }
        this.initialPileup = noOverlapPileup.getPileupWithoutDeletions();
        this.qualityScoreFilteredPileup = initialPileup.getBaseFilteredPileup(minQualityScore);
        this.finalPileup = qualityScoreFilteredPileup.getPileupWithoutMappingQualityZeroReads();

        // apply the same process to derive the strand-specific pileups but in this case don't disregard the overlaps
        this.finalPileupPositiveStrand = pileup.getPileupWithoutDeletions().getBaseFilteredPileup(minQualityScore).getPileupWithoutMappingQualityZeroReads().getPositiveStrandPileup();
        this.finalPileupNegativeStrand = pileup.getPileupWithoutDeletions().getBaseFilteredPileup(minQualityScore).getPileupWithoutMappingQualityZeroReads().getNegativeStrandPileup();

        for (PileupElement p : qualityScoreFilteredPileup) {
            if (p.getMappingQual() == 0 && !allowMapq0ForQualSum) { continue; }
            if (p.getQual() <= minQSumQualityScore) { continue; }

            if (p.getQual() > this.minQSumQualityScore) { qualitySums.incrementSum((char)p.getBase(), 1, p.getQual()); }
 
        }

        this.finalPileupReads = finalPileup.getReads();
        this.finalPileupOffsets = finalPileup.getOffsets();

        // calculate how many are at this site and how many insertions are within INSERTION_PROXIMITY bp
        for (PileupElement p : noOverlapPileup) {
            if (p.getBase() == PileupElement.DELETION_BASE) {
                deletionsCount++;
            } else {

                // check for nearby events
                if (p.getRead().getCigarString().contains("I") || p.getRead().getCigarString().contains("D")) {
                    int eventStart = 0;

                    for (CigarElement ce : p.getRead().getCigar().getCigarElements()) {
                        if (ce.getOperator() == CigarOperator.INSERTION && (Math.abs(eventStart - p.getOffset()) <= GAP_EVENT_PROXIMITY)) {
                            insertionsCount++;
                            break;
                        }

                        if (ce.getOperator() == CigarOperator.DELETION && (Math.abs(eventStart - p.getOffset()) <= GAP_EVENT_PROXIMITY)) {
                            deletionsCount++;
                            break;
                        }

                        eventStart += ce.getLength();
                    }
                }
            }
        }

    }

    public static ReadBackedPileup getOverlappingFragmentFilteredPileupButPreferMismatches(ReadBackedPileup rbp, byte ref) {
        return getOverlappingFragmentFilteredPileup(rbp, ref, true);
    }

    public static ReadBackedPileup getOverlappingFragmentFilteredPileup(ReadBackedPileup rbp, byte ref) {
        return getOverlappingFragmentFilteredPileup(rbp, ref, false);
    }

    public static ReadBackedPileup getOverlappingFragmentFilteredPileup(ReadBackedPileup rbp, byte ref, boolean retainMismatches) {
        Map<String,Integer> filteredPileup = new HashMap<String, Integer>();

        PileupElement[] listOfElements = new PileupElement[rbp.getNumberOfElements()];
        boolean[] elementsToKeep = new boolean[rbp.getNumberOfElements()];
        int i=0;
        for ( PileupElement p : rbp ) {
            listOfElements[i] = p;
            String readName = p.getRead().getReadName();

            // if we've never seen this read before, life is good
            if (!filteredPileup.containsKey(readName)) {
                filteredPileup.put(readName, i);
                elementsToKeep[i] = true;
            } else {
                int existingIndex = filteredPileup.get(readName);
                PileupElement existing = listOfElements[existingIndex];

                // if the reads disagree at this position...
                if (existing.getBase() != p.getBase()) {
                    //... and we're not retaining mismatches, throw them both out
                    if (!retainMismatches) {
                        filteredPileup.remove(readName);
                        elementsToKeep[existingIndex] = false;

                    //... and we are retaining mismatches, keep the mismatching one
                    } else {
                        if (p.getBase() != ref) {
                            elementsToKeep[existingIndex] = false;
                            filteredPileup.put(readName, i);
                            elementsToKeep[i] = true;
                        }
                    }
                // Otherwise, keep the element with the higher quality score
                } else {
                    if (existing.getQual() < p.getQual()) {
                        elementsToKeep[existingIndex] = false;
                        filteredPileup.put(readName, i);
                        elementsToKeep[i] = true;
                    }
                }
            }
            i++;
        }

        List<PileupElement> sortedList = new ArrayList<PileupElement>(rbp.getNumberOfElements());
        for(int j=0; j<elementsToKeep.length; j++) {
            Boolean b = elementsToKeep[j];
            if (b != null && b ) { sortedList.add(listOfElements[j]); }
        }

        return new ReadBackedPileupImpl(rbp.getLocation(), sortedList);
    }

    public double estimateAlleleFraction(char ref, char alt) {
        return estimateAlleleFraction(this.finalPileup, ref, alt);
    }

    public static double estimateAlleleFraction(ReadBackedPileup pileup, char ref, char alt) {
        int[] counts = pileup.getBaseCounts();
        double refCount = (double) counts[BaseUtils.simpleBaseToBaseIndex((byte)ref)];
        double altCount = (double) counts[BaseUtils.simpleBaseToBaseIndex((byte)alt)];
        double depth = refCount + altCount;
        return (depth==0)?0:(altCount / depth);
    }

    public double calculateLogLikelihood(byte alt, double f) {
        return LocusReadPile.calculateLogLikelihood(this.finalPileup, ((byte) this.refBase), alt, f);
    }

    public double calculateAltVsRefLOD(ReadBackedPileup pileup, byte alt, double fAlternate, double fReference) {
        return LocusReadPile.calculateAltVsRefLOD(pileup, ((byte) this.refBase), alt, fAlternate, fReference);
    }

    public static double calculateAltVsRefLOD(ReadBackedPileup pileup, byte ref, byte alt, double fAlternate, double fReference) {
        double lodAlt = LocusReadPile.calculateLogLikelihood(pileup, ref, alt, fAlternate);
        double lodRef = LocusReadPile.calculateLogLikelihood(pileup, ref, alt, fReference);
        return lodAlt - lodRef;
    }

    public double calculateAltVsRefLOD(byte alt, double fAlternate, double fReference) {
        return calculateAltVsRefLOD(this.finalPileup, alt, fAlternate, fReference);
    }

    public double calculateRefVsAltLOD(ReadBackedPileup pileup, byte alt, double fAlternate, double fReference) {
        return -1*calculateAltVsRefLOD(pileup, alt, fAlternate, fReference);
    }

    static public double calculateLogLikelihood(ReadBackedPileup pileup, byte ref, byte alt, double f) {

        double ll = 0;
        for(PileupElement pe : pileup) {
            byte base = pe.getBase();
            byte qual = pe.getQual();
            double e = pow(10, (qual / -10.0));

            if (base == ref) {
                ll += Math.log10(f*e/3 + (1-f)*(1-e));
            } else if (base == alt) {
                ll += Math.log10(f*(1-e) + (1-f)*e/3);
            } else {
                ll += Math.log10(2*e/3);
            }

        }
        return ll;
    }

    public VariableAllelicRatioGenotypeLikelihoods calculateLikelihoods(ReadBackedPileup pileup) {
        return calculateLikelihoods(0.5f, pileup);
    }

    public VariableAllelicRatioGenotypeLikelihoods calculateLikelihoods(double alpha, ReadBackedPileup pileup) {
        // TODO: see if we can move to use "likelihoods.add(pileup, true, true, this.minQualityScore)"
        VariableAllelicRatioGenotypeLikelihoods likelihoods
                = new VariableAllelicRatioGenotypeLikelihoods(refBase, alpha);

        // we have to do this rather than pass in the ReadBackedPileup because that call
        // attempts to make Fragments out of these, which doesn't work if you have
        // a single pileup with multiple samples (as we do in the simulation)
        for(PileupElement pe : pileup) {
            likelihoods.add(pe, false, false, this.minQualityScore);
        }
        return likelihoods;
    }

    public int getFilteredBaseCount(int minBaseQualityScore) {
        return this.finalPileup.getBaseFilteredPileup(minBaseQualityScore).depthOfCoverage();
    }

    public DiploidGenotype getBestGenotype(VariableAllelicRatioGenotypeLikelihoods likelihoods) {
        DiploidGenotype best = null;
        Double bestLikelihood = null;

        for(DiploidGenotype gt : DiploidGenotype.values()) {
            double likelihood = likelihoods.getLikelihood(gt);
            if (bestLikelihood == null || likelihood >= bestLikelihood) {
                best = gt;
                bestLikelihood = likelihood;
            }
        }
        return best;
    }


    public double getHetVsRef(VariableAllelicRatioGenotypeLikelihoods likelihoods, char ref, char altAllele) {
        double[] refHetHom = extractRefHetHom(likelihoods, ref, altAllele);
        return refHetHom[1] - refHetHom[0];
    }

    public double getAltVsRef(VariableAllelicRatioGenotypeLikelihoods likelihoods, char ref, char altAllele) {
        double[] refHetHom = extractRefHetHom(likelihoods, ref, altAllele);
        double tumorLod = logAddSafe(refHetHom[1], refHetHom[2]) - refHetHom[0];

        return tumorLod;
    }

    public static double getRefVsAlt(VariableAllelicRatioGenotypeLikelihoods likelihoods, char ref, char altAllele) {
        double[] refHetHom = extractRefHetHom(likelihoods, ref, altAllele);
        double normalLod = refHetHom[0] - logAddSafe(refHetHom[1], refHetHom[2]);

        return normalLod;
    }


    // ------------------------------------------------------------
    protected static double[] extractRefHetHom(VariableAllelicRatioGenotypeLikelihoods gl, char refAllele, char altAllele) {
        double ref = 0;
        double het = 0;
        double hom = 0;

        for(DiploidGenotype gt : DiploidGenotype.values()) {
            double likelihood = gl.getLikelihood(gt);

            // the reference theory (ref/ref)
            if ( gt.base1 == refAllele && gt.base2 == refAllele) {
                ref = likelihood;
            }

            // the het (ref/alt or alt/ref) theory
            if ( (gt.base1 == refAllele && gt.base2 == altAllele) ||
                 (gt.base1 == altAllele && gt.base2 == refAllele) ) {
                het = likelihood;
            }

            // the hom theory (alt/alt)
            if ( gt.base1 == altAllele && gt.base2 == altAllele) {
                hom = likelihood;
            }
        }

        return new double[]{ref, het, hom};
    }

    /**
     * underflow safe version to compute log(10^a + 10^b)
     *
     *    e.g. compute c as in 10^c = 10^a + 10^b where a > b
     *    which is c = log(10^a + 10^b)
     *             c = log(10^a + 10^(b+a-a))
     *             c = log(10^a + 10^a*10^(b-a))
     *             c = log(10^a *(1 + 10^(b-a))
     *             c = log(10^a) + log(1+10^(b-a))
     *             c = a + log(1+10^(b-a))
     *
     * and in this case if b-a is a very small number, the last term is 0 but we've done this without underflowing on the significant term
     * 
     *
     * e.g. a and b are two log-likelihoods to add
     *
     */
    private static double logAddSafe(double inOne, double inTwo) {
        double a = (inOne > inTwo)?inOne:inTwo;
        double b = (inOne > inTwo)?inTwo:inOne;

        return a + Math.log(1 + Math.pow(10, b-a));
    }


    public int getDeletionsCount() {
        return deletionsCount;
    }

    public int getInsertionsCount() {
        return insertionsCount;
    }
}
