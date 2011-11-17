package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;

import java.util.HashMap;

public class PValuePowerCalculator {
    private static double Q30 = Math.pow(10, (-30/10));
    private static double Q20 = Math.pow(10, (-20/10));

    public static void main(String[] argv) throws MathException {
        double wiggle = 0.000001;
        double fdr = Math.pow(10, -6);

//        Assert.assertEquals(testCalculatePower(50, Q30, fdr, 0.15), 0.955279 , wiggle, "Failed!");
//        Assert.assertEquals(testCalculatePower(20, Q30, fdr, 0.15), 0.5679002 , wiggle, "Failed!");
//        Assert.assertEquals(testCalculatePower(10, Q30, fdr, 0.15), 0.1852461 , wiggle, "Failed!");
//        Assert.assertEquals(testCalculatePower(50, Q30, fdr, 0.35), 0.9999985 , wiggle, "Failed!");
//        Assert.assertEquals(testCalculatePower(50, Q30, fdr, 0.05), 0.2487401 , wiggle, "Failed!");
//        Assert.assertEquals(testCalculatePower(50, Q20, fdr, 0.05), 0.01258744 , wiggle, "Failed!");
//        Assert.assertEquals(testCalculatePower(50, Q20, Math.pow(10, -5), 0.05), 0.03549362 , wiggle, "Failed!");
//
//        // make an instance to test caching
//        PValuePowerCalculator pc = new PValuePowerCalculator(Q30, fdr);
//        Assert.assertEquals(pc.cachingPowerCalculation(50, 0.15), 0.955279 , wiggle, "Failed!");
//        Assert.assertEquals(pc.cachingPowerCalculation(50, 0.15), 0.955279 , wiggle, "Failed!");
//        Assert.assertEquals(pc.cachingPowerCalculation(50, 0.15), 0.955279 , wiggle, "Failed!");
        
    }


    
    private static double testCalculatePower(int n, double eps, double fdr, double delta) throws MathException {
        double power = calculatePower(n, eps, fdr, delta);
        System.out.println("Depth: " + n + " EPS: " + eps + " FDR: " + fdr + " Delta: " + delta + " --> Power: " + power);
        return power;
    }

    private static void log(String s) {
        System.out.println(s);
    }


    private double constantEps;
    private double constantFdr;

    public PValuePowerCalculator(double constantEps, double constantFdr) {
        this.constantEps = constantEps;
        this.constantFdr = constantFdr;
    }

    private static class PowerCacheKey {
        private int n;
        private double delta;

        private PowerCacheKey(int n, double delta) {
            this.n = n;
            this.delta = delta;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            PowerCacheKey that = (PowerCacheKey) o;

            if (Double.compare(that.delta, delta) != 0) return false;
            if (n != that.n) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result;
            long temp;
            result = n;
            temp = delta != +0.0d ? Double.doubleToLongBits(delta) : 0L;
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            return result;
        }
    }

    private HashMap<PowerCacheKey, Double> cache = new HashMap<PowerCacheKey, Double>();

    private double cachingPowerCalculation(int n, double delta) throws MathException {
        PowerCacheKey key = new PowerCacheKey(n, delta);
        Double power = cache.get(key);
        if (power == null) {
            power = calculatePower(n, constantEps, constantFdr, delta);
            cache.put(key, power);
        } else {
            log("Cache hit!");
        }
        return power;        
    }

    private static double calculatePower(int n, double eps, double fdr, double delta) throws MathException {
        BinomialDistribution bdEps = new BinomialDistributionImpl(n, eps);

        double[] p = new double[n+2];
        double[] pv = new double[n+2];

        p[0] = 0;
        for (int i=1; i<=n+1; i++ ) {
            p[i] = bdEps.probability(i-1);

            //log("P("+i+") == " + p[i]);

            // perform the cumulative sum
            pv[i] = pv[i-1] + p[i];
        }

        // log("-------- Step 2 ----------");

        // invert it
        for (int i=0; i<pv.length; i++) {
            pv[i] = 1d - pv[i];
            //log("PV("+i+") == " + pv[i]);
        }

        // log("-------- Step 3 ----------");
        int ks = -1;
        for (int i=0; i<pv.length; i++) {
            if (pv[i] <= fdr) {
                ks = i;
                break;
            }
        }

        if (ks == -1) { throw new MathException("could not compute ks!"); }
        //log("calculated ks == " + ks);

        //log("-------- Step 4 ----------");
        double cval = (fdr - pv[ks] ) / ( pv[ks-1] - pv[ks] );
        //log("calculated cval == " + cval);


        //log("-------- Step 5 ----------");
        BinomialDistribution bdDelta = new BinomialDistributionImpl(n, delta);

        double[] p1 = new double[n+1];

        for (int i=0; i<=n; i++ ) {
            p1[i] = bdDelta.probability(i);

            //log("P1("+i+") == " + p1[i]);
        }

        //log("-------- Step 6 ----------");
        double p1Sum = 0;
        for (int i=0; i<=ks-1; i++) {
            p1Sum += p1[i];
        }
        //log("calculated p1Sum == " + p1Sum);

        double pow = 1 - p1Sum + cval * p1[ks-1];
        //log("calculated power == " + pow);



        return pow;
    }
}
