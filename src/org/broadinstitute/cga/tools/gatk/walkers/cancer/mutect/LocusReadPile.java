package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

import java.util.*;

import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotype;
import org.broadinstitute.sting.utils.pileup.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import static java.lang.Math.pow;

public class LocusReadPile {
    protected ReadBackedPileup pileup;
    List<GATKSAMRecord> finalPileupReads;
    List<Integer> finalPileupOffsets;

    protected char refBase;
    protected int minQualityScore;
    protected int minQSumQualityScore;
    protected QualitySums qualitySums = new QualitySums();

    public ReadBackedPileup initialPileup;
    protected ReadBackedPileup qualityScoreFilteredPileup;
    public ReadBackedPileup finalPileup;


    public ReadBackedPileup finalPileupPositiveStrand;
    public ReadBackedPileup finalPileupNegativeStrand;

    public static final int GAP_EVENT_PROXIMITY = 5; // a 11bp window
    protected int deletionsCount = 0;
    protected int insertionsCount = 0;

    public LocusReadPile(ReadBackedPileup pileup, char refBase, int minQualityScore, int minQSumQualityScore) {
        this(pileup, refBase, minQualityScore, minQSumQualityScore, false, false);
    }

    public LocusReadPile(ReadBackedPileup pileup, char refBase, int minQualityScore, int minQSumQualityScore, boolean allowMapq0ForQualSum, boolean retainOverlapMismatches) {
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

            if (p.getQual() > this.minQSumQualityScore) { qualitySums.incrementSum((char)p.getBase(), p.getRepresentativeCount(), p.getRepresentativeCount() * p.getQual()); }
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
        double lodAlt = LocusReadPile.calculateLogLikelihood(pileup, ((byte) this.refBase), alt, fAlternate);
        double lodRef = LocusReadPile.calculateLogLikelihood(pileup, ((byte) this.refBase), alt, fReference);
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
     * @return
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
