package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

import java.util.*;

import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotype;
import org.broadinstitute.sting.utils.pileup.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

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

    public static final int INSERTION_EVENT_PROXIMITY = 5; // a 11bp window
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


        for (PileupElement p : qualityScoreFilteredPileup) {
            if (p.getMappingQual() == 0 && !allowMapq0ForQualSum) { continue; }
            if (p.getQual() <= minQSumQualityScore) { continue; }

            if (p.getQual() > this.minQSumQualityScore) { qualitySums.incrementSum((char)p.getBase(), p.getQual()); }
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
                        if (ce.getOperator() == CigarOperator.INSERTION && (Math.abs(eventStart - p.getOffset()) <= INSERTION_EVENT_PROXIMITY)) {
                            insertionsCount++;
                            break;
                        }

                        if (ce.getOperator() == CigarOperator.DELETION && (Math.abs(eventStart - p.getOffset()) <= INSERTION_EVENT_PROXIMITY)) {
                            deletionsCount++;
                            break;
                        }

                        eventStart += ce.getLength();
                    }
                }
            }
        }

        // get the number of deletions
//        this.deletionsCount = noOverlapPileup.getNumberOfDeletions();
        
    }

    public static ReadBackedPileup getOverlappingFragmentFilteredPileupButPreferMismatches(ReadBackedPileup rbp, byte ref) {
        return getOverlappingFragmentFilteredPileup(rbp, ref, true);
    }

    public static ReadBackedPileup getOverlappingFragmentFilteredPileup(ReadBackedPileup rbp, byte ref) {
        return getOverlappingFragmentFilteredPileup(rbp, ref, false);
    }

    public static ReadBackedPileup getOverlappingFragmentFilteredPileup(ReadBackedPileup rbp, byte ref, boolean retainMismatches) {
        Map<String,PileupElement> filteredPileup = new HashMap<String, PileupElement>();

        for ( PileupElement p : rbp ) {
            String readName = p.getRead().getReadName();

            // if we've never seen this read before, life is good
            if (!filteredPileup.containsKey(readName)) {
                filteredPileup.put(readName, p);
            } else {
                PileupElement existing = filteredPileup.get(readName);

                // if the reads disagree at this position...
                if (existing.getBase() != p.getBase()) {
                    //... and we're not retaining mismatches, throw them both out
                    if (!retainMismatches) {
                        filteredPileup.remove(readName);

                    //... and we are retaining mismatches, keep the mismatching one
                    } else {
                        if (p.getBase() != ref) {
                            filteredPileup.put(readName, p);
                        }
                    }
                // Otherwise, keep the element with the higher quality score
                } else {
                    if (existing.getQual() < p.getQual()) {
                        filteredPileup.put(readName, p);
                    }
                }
            }
        }

        // TODO: hacked so that we returned a RBPL constructed with a sorted list of elements.  should really be using the native fragment stuff throughout...
        SortedSet<PileupElement> filteredPileupElements = new TreeSet<PileupElement>(new PileupElementAlignmentStartPositionComparator());
        for(PileupElement filteredElement: filteredPileup.values()) {
            filteredPileupElements.add(filteredElement);
        }

        List<PileupElement> sortedList = new ArrayList<PileupElement>(filteredPileupElements.size());
        for(PileupElement pe : filteredPileupElements) {
            sortedList.add(pe);
        }

        return new ReadBackedPileupImpl(rbp.getLocation(), sortedList);
    }

    private static class PileupElementAlignmentStartPositionComparator implements Comparator<PileupElement>{
        public int compare(PileupElement pe1, PileupElement pe2) {
            if (pe1.getRead().getAlignmentStart() < pe2.getRead().getAlignmentStart()) { return -1; } else { return 1; }
        }
    }

    public double estimateAlleleFraction(char ref, char alt) {
        return estimateAlleleFraction(this.finalPileup, ref, alt);
    }

    public static double estimateAlleleFraction(ReadBackedPileup pileup, char ref, char alt) {
        int[] counts = pileup.getBaseCounts();
        double refCount = (double) counts[BaseUtils.simpleBaseToBaseIndex(ref)];
        double altCount = (double) counts[BaseUtils.simpleBaseToBaseIndex(alt)];
        double depth = refCount + altCount;
        return (depth==0)?0:(altCount / depth);
    }

    public VariableAllelicRatioGenotypeLikelihoods calculateLikelihoods(ReadBackedPileup pileup) {
        return calculateLikelihoods(0.5f, pileup);
    }

    public VariableAllelicRatioGenotypeLikelihoods calculateLikelihoods(double alpha, ReadBackedPileup pileup) {
        VariableAllelicRatioGenotypeLikelihoods likelihoods
                = new VariableAllelicRatioGenotypeLikelihoods(refBase, alpha);
        likelihoods.add(pileup, false, false, this.minQualityScore);
        return likelihoods;
    }

    public int getFilteredBaseCount(int minBaseQualityScore) {
        return this.finalPileup.getBaseFilteredPileup(minBaseQualityScore).depthOfCoverage();
    }

    public List<Byte> getLocusBases(int locusOffset) {
        List<Byte> bases = new ArrayList<Byte>(finalPileupReads.size());

        for(int i=0; i< finalPileupReads.size(); i++) {
            SAMRecord read = finalPileupReads.get(i);
            int readOffset = finalPileupOffsets.get(i);

            int offset = readOffset + locusOffset;
            if (offset >= 0 && offset < read.getReadString().length()) {
                byte base = read.getReadBases()[offset];
                if (base != (byte) 'N') {
                    bases.add(read.getReadBases()[offset]);
                }
            }
        }
        return bases;
    }

//    public List<Byte> getLocusQuals(int locusOffset) {
//        List<Byte> quals = new ArrayList<Byte>(finalPileup.getReads().size());
//
//        for(int i=0; i<finalPileup.getReads().size(); i++) {
//            SAMRecord read = finalPileup.getReads().get(i);
//            int readOffset = finalPileup.getOffsets().get(i);
//
//            int offset = readOffset + locusOffset;
//            if (offset >= 0 && offset < read.getReadString().length()) {
//                quals.add(read.getBaseQualities()[offset]);
//            }
//        }
//        return quals;
//    }
//

    public List<Byte> getQualityScores(char allele) {
        return getQualityScores(0, allele);
    }

    public List<Byte> getQualityScores(int locusOffset, char allele) {
        List<Byte> scores = new ArrayList<Byte>();

        for(int i=0; i<finalPileupReads.size(); i++) {
            SAMRecord read = finalPileupReads.get(i);
            int readOffset = finalPileupOffsets.get(i);

            int offset = readOffset + locusOffset;
            if (offset >= 0 && offset < read.getReadString().length()) {
                char base = read.getReadString().charAt(offset);
                byte qual = read.getBaseQualities()[offset];

                if (base == allele) { scores.add(qual); }
            }
        }
        return scores;        
    }

    public List<Integer> getPositionsInRead(char allele) {
        return getPositionsInRead(0, allele);
    }

    public List<Integer> getPositionsInRead(int locusOffset, char allele) {
        List<Integer> positions = new ArrayList<Integer>();

        for(int i=0; i<finalPileupReads.size(); i++) {
            SAMRecord read = finalPileupReads.get(i);
            int readOffset = finalPileupOffsets.get(i);

            int offset = readOffset + locusOffset;
            if (offset >= 0 && offset < read.getReadString().length()) {
                char base = read.getReadString().charAt(offset);
                if (base == allele) { positions.add(offset); }
            }
        }
        return positions;
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

    public static double getRefVsNextBest(VariableAllelicRatioGenotypeLikelihoods likelihoods, char ref) {
        Double refLikelihood = null;
        Double nextBest = null;

        for(DiploidGenotype gt : DiploidGenotype.values()) {
            double likelihood = likelihoods.getLikelihood(gt);

            // the reference theory (ref/ref)
            if ( gt.base1 == ref && gt.base2 == ref) {
                refLikelihood = likelihood;
            } else {
                if (nextBest == null || likelihood >= nextBest) {
                    nextBest = likelihood; 
                }
            }

        }

        return refLikelihood - nextBest;
    }

    public double getHetVsRef(VariableAllelicRatioGenotypeLikelihoods likelihoods, char ref, char altAllele) {
        double[] refHetHom = extractRefHetHom(likelihoods, ref, altAllele);
        return refHetHom[1] - refHetHom[0];
    }

    public double getRefVsHet(VariableAllelicRatioGenotypeLikelihoods likelihoods, char ref, char altAllele) {
        double[] refHetHom = extractRefHetHom(likelihoods, ref, altAllele);
        return refHetHom[0] - refHetHom[1];
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

    /**
     * Extract the LOD comparing ref:ref to ref:alt and alt:alt
     */
    private double[] extractRefAlt(VariableAllelicRatioGenotypeLikelihoods gl, char ref, char altAllele) {
        double refRef = 0;
        double altRef = 0;
        double altAlt = 0;

        // FIXME: potential underflow problem!        
        for(DiploidGenotype gt : DiploidGenotype.values()) {
            double likelihood = gl.getLikelihood(gt);

            // the ref:mutant theory
            if ( (gt.base1 == ref && gt.base2 == altAllele) ||
                 (gt.base1 == altAllele && gt.base2 == ref) ) {
                altRef += Math.pow(10, likelihood);
            }

            if ( gt.base1 == altAllele && gt.base2 == altAllele) {
                altAlt += Math.pow(10, likelihood);
            }

            if ( gt.base1 == ref && gt.base2 == ref) {
                refRef = likelihood;
            }

        }
        return new double[]{refRef, altRef, altAlt};
    }

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
