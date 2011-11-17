package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;

public class NormalPowerCalculator extends AbstractPowerCalculator {

    private static double Q30_EPS = Math.pow(10, (-30/10)) / 3d;
    private static double LOD = 2.3d;

    public static void main(String[] argv) throws MathException {
        double wiggle = 0.000001;

        test(testCalculatePower(50, Q30_EPS, LOD), 1.0000000, wiggle, "Failed!");
        test(testCalculatePower(10, Q30_EPS, LOD), 0.9973492, wiggle, "Failed!");
        test(testCalculatePower( 8, Q30_EPS, LOD), 0.9974184, wiggle, "Failed!");
        test(testCalculatePower( 5, Q30_EPS, LOD), 0.0000000, wiggle, "Failed!");
        test(testCalculatePower(10, Math.pow(10, (-20/10)) / 3d, LOD), 0.9762534, wiggle, "Failed!");


        // make an instance to test caching
        NormalPowerCalculator pc = new NormalPowerCalculator(Q30_EPS, LOD);
        test(pc.cachingPowerCalculation(10), 0.9973492, wiggle, "Failed!");
        test(pc.cachingPowerCalculation(10), 0.9973492, wiggle, "Failed!");
        test(pc.cachingPowerCalculation(10), 0.9973492, wiggle, "Failed!");

    }

    private static double testCalculatePower(int n, double eps, double threshold) throws MathException {
        double power = calculatePower(n, eps, threshold);
        System.out.println("Depth: " + n + " EPS: " + eps + " Threshols: " + threshold + " --> Power: " + power);
        return power;
    }

    public NormalPowerCalculator(double constantEps, double constantLodThreshold) {
        this.constantEps = constantEps;
        this.constantLodThreshold = constantLodThreshold;
    }



    public double cachingPowerCalculation(int n) throws MathException {
        PowerCacheKey key = new PowerCacheKey(n, 0.5);
        Double power = cache.get(key);
        if (power == null) {
            power = calculatePower(n, constantEps, constantLodThreshold);
            cache.put(key, power);
        }
        return power;
    }

    private static double calculateNormalLod(int depth, int alts, double eps) {
        return (calculateLogLikelihood(depth, alts, eps, 0) - calculateLogLikelihood(depth, alts, eps, 0.5));
    }

    private static double calculatePower(int depth, double eps, double lodThreshold) throws MathException {
        if (depth==0) return 0;

        // calculate the probability of each configuration
		// NOTE: reorder from tumor case to be # of ref instead of # of alts
        BinomialDistribution binom = new BinomialDistributionImpl(depth, eps);
        double[] p = new double[depth+1];
        for(int i=0; i<p.length; i++) {
            p[i] = binom.probability(depth - i);
        }

        // calculate the LOD scores
        double[] lod = new double[depth+1];
        for(int i=0; i<lod.length; i++) {
            lod[i] = calculateNormalLod(depth, depth-i, eps);
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

        double x = 1d - (lodThreshold - lod[k-1]) / (lod[k] - lod[k-1]);
        double power = x*p[k-1];
        for(int i=k; i<p.length; i++) {
            power += p[i];
        }

        return(power);
    }

}
