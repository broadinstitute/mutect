package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import org.broadinstitute.sting.utils.BaseUtils;

import static java.lang.Math.log10;
import static java.lang.Math.pow;

public class VariableAllelicRatioGenotypeLikelihoods extends DiploidSNPGenotypeLikelihoods {
    protected double logF;
    protected double logOneMinusF;
    protected double logHalf;
    protected char ref;


    /**
     * Create a new GenotypeLikelhoods object with flat priors for each diploid genotype
     *
     * @param ref reference base
     * @param f non-reference allele fraction estimate
     */
    public VariableAllelicRatioGenotypeLikelihoods(char ref, double f) {
//        super(new DiploidSNPGenotypePriors(), DiploidSNPGenotypeLikelihoods.DEFAULT_PCR_ERROR_RATE);
        // TODO: re-enable the non-zero PCR Error rate once we understand the effect of it
        super(new DiploidSNPGenotypePriors(), 0);
        this.ref = ref;

        this.logF = log10(1/f);
        this.logOneMinusF = log10(1/(1-f));
        this.logHalf = log10(1/.5);
    }

    protected DiploidSNPGenotypeLikelihoods calculateGenotypeLikelihoods(byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2) {
        double[] log10FourBaseLikelihoods = computeLog10Likelihoods(observedBase1, qualityScore1, observedBase2, qualityScore2);

        try {

            VariableAllelicRatioGenotypeLikelihoods gl = (VariableAllelicRatioGenotypeLikelihoods)this.clone();
            gl.setToZero();

            for ( DiploidGenotype g : DiploidGenotype.values() ) {

                double fBase1;
                double fBase2;
                if (g.base1 == ref || g.base2 == ref) {

                    // if it is ref/ref use half
                    if (g.base1 == g.base2) {
                        fBase1 = logHalf;
                        fBase2 = logHalf;
                    // if one base is reference, use f
                    } else {
                        fBase1 = (g.base1 == ref)? logOneMinusF : logF;
                        fBase2 = (g.base2 == ref)? logOneMinusF : logF;
                    }
                } else {
                    fBase1 = logHalf;
                    fBase2 = logHalf;
                }

                double p_base = 0.0;
                p_base += pow(10, log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(g.base1)] - fBase1);
                p_base += pow(10, log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(g.base2)] - fBase2);
                double likelihood = log10(p_base);

                gl.log10Likelihoods[g.ordinal()] += likelihood;
                gl.log10Posteriors[g.ordinal()] += likelihood;
            }
            if ( VERBOSE ) {
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%s\t", g); }
                System.out.println();
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", gl.log10Likelihoods[g.ordinal()]); }
                System.out.println();
            }

            return gl;

         } catch ( CloneNotSupportedException e ) {
             throw new RuntimeException(e);
         }
    }

    public int add(byte observedBase1, byte qualityScore1) {
        byte observedBase2 = 0, qualityScore2 = 0;

        // Just look up the cached result if it's available, or compute and store it
        DiploidSNPGenotypeLikelihoods gl;
        if ( ! inCache(observedBase1, qualityScore1, observedBase2, qualityScore2, FIXED_PLOIDY) ) {
            gl = calculateCachedGenotypeLikelihoods(observedBase1, qualityScore1, observedBase2, qualityScore2, FIXED_PLOIDY);
        } else {
            gl = getCachedGenotypeLikelihoods(observedBase1, qualityScore1, observedBase2, qualityScore2, FIXED_PLOIDY);
        }

        double[] likelihoods = gl.getLikelihoods();

        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double likelihood = likelihoods[g.ordinal()];

            log10Likelihoods[g.ordinal()] += likelihood;
            log10Posteriors[g.ordinal()] += likelihood;
        }

        return 1;
    }


    ////////////////////////////////////////////////////////////////////////////////////////
    //
    // Need to disable caching for now since it's not built to handle variable allelic fractions
    //
    ////////////////////////////////////////////////////////////////////////////////////////

    @Override
    protected boolean inCache(byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2, int ploidy) {
        return false;
    }

    @Override
    protected void setCache( DiploidSNPGenotypeLikelihoods[][][][][] cache,
                             byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2, int ploidy,
                             DiploidSNPGenotypeLikelihoods val ) {
        // do nothing
    }
}
