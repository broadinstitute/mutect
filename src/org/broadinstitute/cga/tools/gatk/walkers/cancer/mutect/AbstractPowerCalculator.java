package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import java.util.HashMap;

public class AbstractPowerCalculator {
    protected HashMap<PowerCacheKey, Double> cache = new HashMap<PowerCacheKey, Double>();
    protected double constantEps;
    protected double constantLodThreshold;

    protected static class PowerCacheKey {
        private int n;
        private double delta;

        public PowerCacheKey(int n, double delta) {
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

    protected static double calculateLogLikelihood(int depth, int alts, double eps, double f) {
        double a = (depth-alts) * Math.log10(f*eps + (1d-f)*(1d-eps));
        double b = (alts) * Math.log10(f*(1d-eps) + (1d-f)*eps);
        return (a+b);
    }

    protected static void test(double actual, double expected, double tolerance, String msg) {
        if (Math.abs(actual - expected) > tolerance) {
            throw new RuntimeException("FAILED: calculated " + actual + " but expected " + expected + " -> " + msg);
        }
    }

}
