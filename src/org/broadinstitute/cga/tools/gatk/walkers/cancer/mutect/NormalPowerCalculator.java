package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;

public class NormalPowerCalculator extends AbstractPowerCalculator {
    private boolean constantEnableSmoothing;

    public NormalPowerCalculator(double constantEps, double constantLodThreshold) {
        this.constantEps = constantEps;
        this.constantLodThreshold = constantLodThreshold;
        this.constantEnableSmoothing = true;
    }

    public NormalPowerCalculator(double constantEps, double constantLodThreshold, boolean enableSmoothing) {
        this.constantEps = constantEps;
        this.constantLodThreshold = constantLodThreshold;
        this.constantEnableSmoothing = enableSmoothing;
    }

    public double cachingPowerCalculation(int n) throws MathException {
        PowerCacheKey key = new PowerCacheKey(n, 0.5);
        Double power = cache.get(key);
        if (power == null) {
            power = calculatePower(n, constantEps, constantLodThreshold, constantEnableSmoothing);
            cache.put(key, power);
        }
        return power;
    }

    protected static double calculateNormalLod(int depth, int alts, double eps) {
        return (calculateLogLikelihood(depth, alts, eps, 0) - calculateLogLikelihood(depth, alts, eps, 0.5));
    }

    protected static double calculatePower(int depth, double eps, double lodThreshold, boolean enableSmoothing) throws MathException {
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

        double power = 0;

        // here we correct for the fact that the exact lod threshold is likely somewhere between
        // the k and k-1 bin, so we prorate the power from that bin
        // the k and k-1 bin, so we prorate the power from that bin
        // if k==0, it must be that lodThreshold == lod[k] so we don't have to make this correction
        if ( enableSmoothing && k > 0 ) {
            double x = 1d - (lodThreshold - lod[k-1]) / (lod[k] - lod[k-1]);
            power = x*p[k-1];
        }

        for(int i=k; i<p.length; i++) {
            power += p[i];
        }

        return(power);
    }

}
