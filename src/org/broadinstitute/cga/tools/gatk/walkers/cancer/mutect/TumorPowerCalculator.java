package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;

import java.util.HashMap;

public class TumorPowerCalculator extends AbstractPowerCalculator{

    private static double Q30_EPS = Math.pow(10, (-30/10)) / 3d;
    private static double LOD = 6.3d;

    public static void main(String[] argv) throws MathException {
        double wiggle = 0.00001;

        // test log likelihood calculation
        test(calculateLogLikelihood(50, 5, 0.001, 0.1), -7.059164, wiggle, "Failed!");
        test(calculateLogLikelihood(40, 5, 0.001, 0.1), -6.597727, wiggle, "Failed!");
        test(calculateLogLikelihood(40, 10, 0.001, 0.1), -11.34971, wiggle, "Failed!");
        test(calculateLogLikelihood(40, 10, 0.01, 0.1), -11.15482, wiggle, "Failed!");
        test(calculateLogLikelihood(40, 10, 0.01, 0.2), -9.866713, wiggle, "Failed!");


        test(testCalculatePower(50, Q30_EPS, LOD, 0.15, 0), 0.975342, wiggle, "Failed!");

        test(testCalculatePower(20, Q30_EPS, LOD, 0.15, 0), 0.6364833, wiggle, "Failed!");
        test(testCalculatePower(10, Q30_EPS, LOD, 0.15, 0), 0.3167154, wiggle, "Failed!");
        test(testCalculatePower(50, Q30_EPS, LOD, 0.35, 0), 0.9999994, wiggle, "Failed!");
        test(testCalculatePower(50, Q30_EPS, LOD, 0.05, 0), 0.3893047, wiggle, "Failed!");
        test(testCalculatePower(50, Q30_EPS, 4.3, 0.05, 0), 0.6068266, wiggle, "Failed!");
        test(testCalculatePower(50, Q30_EPS, 4.3, 0.05, 0.02), 0.0039135, wiggle, "Failed!");

        // make an instance to test caching
        TumorPowerCalculator pc = new TumorPowerCalculator(Q30_EPS, LOD, 0);
        test(pc.cachingPowerCalculation(50, 0.15), 0.975342, wiggle, "Failed!");
        test(pc.cachingPowerCalculation(50, 0.15), 0.975342, wiggle, "Failed!");
        test(pc.cachingPowerCalculation(50, 0.15), 0.975342, wiggle, "Failed!");

    }

    private static double testCalculatePower(int n, double eps, double fdr, double delta, double contam) throws MathException {
        double power = calculatePower(n, eps, fdr, delta, contam);
        System.out.println("Depth: " + n + " EPS: " + eps + " FDR: " + fdr + " Delta: " + delta + " --> Power: " + power);
        return power;
    }


    private double constantContamination;

    public TumorPowerCalculator(double constantEps, double constantLodThreshold, double constantContamination) {
        this.constantEps = constantEps;
        this.constantLodThreshold = constantLodThreshold;
        this.constantContamination = constantContamination;
    }


    public double cachingPowerCalculation(int n, double delta) throws MathException {
        PowerCacheKey key = new PowerCacheKey(n, delta);
        Double power = cache.get(key);
        if (power == null) {
            power = calculatePower(n, constantEps, constantLodThreshold, delta, constantContamination);
            cache.put(key, power);
        }
        return power;        
    }




    private static double calculateTumorLod(int depth, int alts, double eps, double contam) {
        double f = (double) alts / (double) depth;
        return (calculateLogLikelihood(depth, alts, eps, f) - calculateLogLikelihood(depth, alts, eps, Math.min(f,contam)));
    }

    private static double calculatePower(int depth, double eps, double lodThreshold, double delta, double contam) throws MathException {
        if (depth==0) return 0;

        int n = depth; // for consistency with R code which has an extra outer loop


        // calculate the probability of each configuration
        double p_alt_given_e_delta = delta*(1d-eps) + (1d-delta)*eps;
        BinomialDistribution binom = new BinomialDistributionImpl(n, p_alt_given_e_delta);
        double[] p = new double[depth+1];
        for(int i=0; i<p.length; i++) {
            p[i] = binom.probability(i);
        }

        // calculate the LOD scores
        double[] lod = new double[depth+1];
        for(int i=0; i<lod.length; i++) {
            lod[i] = calculateTumorLod(n, i, eps, contam);
        }

        // TODO: optimization -- either solve directly, or take advantage of the fact that LODs monotonically increase to find the threshold and then just do 1-cumsum
        int k = -1;
        for(int i=0; i<lod.length; i++) {
            if (lod[i] >= lodThreshold) {
                k = i;
                break;
            }
        }

        // if no depth meets the lod score, the power is zero
        if (k == -1) {
            return 0;
        }

        double power = 0;

        // here we correct for the fact that the exact lod threshold is likely somewhere between
        // the k and k-1 bin, so we prorate the power from that bin
        // the k and k-1 bin, so we prorate the power from that bin
        // if k==0, it must be that lodThreshold == lod[k] so we don't have to make this correction
        if ( k > 0 ) {
            double x = 1d - (lodThreshold - lod[k-1]) / (lod[k] - lod[k-1]);
            power = x*p[k-1];
        }

        for(int i=k; i<p.length; i++) {
            power += p[i];
        }

        return(power);
    }

}
