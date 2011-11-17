package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import cern.jet.random.Normal;
import cern.jet.random.AbstractDistribution;

import java.util.*;

import cern.jet.random.engine.MersenneTwister;
import net.sf.picard.util.Histogram;

/**
 * Implementation of a wilcoxon rank sum test
 *
 * @author Tim Fennell
 */
public class RankSumTest {
    /**
     * Simple class to hold a value, it's source series, and it's rank.
     */
    private static final class Rank implements Comparable<Rank> {
        final int value;
        float rank;
        final int series;

        private Rank(int value, float rank, int series) {
            this.value = value;
            this.rank = rank;
            this.series = series;
        }

        public int compareTo(Rank that) {
            return this.value - that.value;
        }

        @Override public String toString() {
            return "Rank{" +
                    "value=" + value +
                    ", rank=" + rank +
                    ", series=" + series +
                    '}';
        }
    }

    /**
     * The results of performing a rank sum test.
     */
    public static class Result {
        private final double u;
        private final double z;
        private final double p;
        private final double medianShift;

        public Result(double u, double z, double p, double medianShift) {
            this.u = u;
            this.z = z;
            this.p = p;
            this.medianShift = medianShift;
        }

        public double getU() { return u; }
        public double getZ() { return z; }
        public double getP() { return p; }
        public double getMedianShift() { return medianShift; }
    }

    // Constructs a normal distribution; actual values of mean and SD don't matter since it's
    // just used to convert a z-score into a cumulative probability
    private static final double NORMAL_MEAN = 100;
    private static final double NORMAL_SD   = 15;
    private static final Normal NORMAL = new Normal(NORMAL_MEAN, NORMAL_SD, null);

    /**
     * The minimum length for both data series (individually) in order to use a normal distribution
     * to calculate Z and p. If either series is shorter than this value then a custom distribution
     * will be used.
     */
    private int minimumNormalN = 20;

    /**
     * When generating custom distributions how many random tests should be performed to
     * construct the distribution.
     */    
    private int valuesInCustomDistribution = 100000;

    /**
     * Distribution from which values are sampled when constructing test data sets to build
     * distributions of U for small sample sizes.
     */
    private static int SEED = 32674;
    private AbstractDistribution modelDistribution = new Normal(22, 6, new MersenneTwister(SEED));


    /**
     * A map of (n1,n2)->distribution that can be used to calculate z and p for values of u
     * with not enough data points to approximate the normal distribution
     */
    private final Map<Key,Histogram<Double>> SMALL_DISTRIBUTIONS = new HashMap<Key,Histogram<Double>>();

    /** Sets the minimum number of values in each data series to use the normal distribution approximation. */
    public void setMinimumSeriesLengthForNormalApproximation(final int n) {
        this.minimumNormalN = n;
    }

    /** Sets how many random tests should be performed when constructing custom distributions. */
    public void setValuesInCustomDistribution(int valuesInCustomDistribution) {
        this.valuesInCustomDistribution = valuesInCustomDistribution;
    }

    /**
     * Sets the model distribution used for creating simulated datasets to build distributions
     * of U for small sample sizes.
     */
    public void setModelDistribution(AbstractDistribution modelDistribution) {
        this.modelDistribution = modelDistribution;
    }

    /**
     * Calculates the rank-sum test statisic U (sometimes W) from two sets of input data
     */
    public double calculateU(final int[] series1, final int[] series2) {
        Arrays.sort(series1);
        Arrays.sort(series2);

        // Make a merged ranks array
        final Rank[] ranks = new Rank[series1.length + series2.length];
        {
            int i=0, j=0, r=0;
            while (r < ranks.length) {
                if (i >= series1.length) {
                    ranks[r++] = new Rank(series2[j++], r, 2);
                }
                else if (j >= series2.length) {
                    ranks[r++] = new Rank(series1[i++], r, 1);
                }
                else if (series1[i] <= series2[j]) {
                    ranks[r++] = new Rank(series1[i++], r, 1);
                }
                else {
                    ranks[r++] = new Rank(series2[j++], r, 2);
                }
            }
        }

        // Now sort out any tie bands
        for (int i=0; i<ranks.length;) {
            float rank = ranks[i].rank;
            int count = 1;
            for (int j=i+1; j<ranks.length && ranks[j].value == ranks[i].value; ++j) {
                rank += ranks[j].rank;
                ++count;
            }

            if (count > 1) {
                rank /= count;
                for (int j=i; j<i+count; ++j) {
                    ranks[j].rank = rank;
                }
            }

            // Skip forward the right number of items
            i += count;
        }

        // Calculate R1 and R2 and U.
        float r1=0, r2=0;
        for (Rank rank : ranks) {
            if (rank.series == 1) r1 += rank.rank;
            else r2 += rank.rank;
        }

        double n1 = series1.length;
        double n2 = series2.length;
        double u1 = r1 - ((n1 * (n1+1)) / 2);
        double u2 = r2 - ((n2 * (n2+1)) / 2);
        return Math.min(u1, u2);
    }

    /**
     * Calculates the Z score (i.e. standard deviations from the mean) of the rank sum
     * test statistics given input data of lengths n1 and n2 respectively.
     */
    public double calculateZ(final double u, final int n1, final int n2) {
        if (n1 >= this.minimumNormalN && n2 >= this.minimumNormalN) {
            double m = (n1 * n2) / 2d;
            double sigma = Math.sqrt((n1*n2*(n1+n2+1)) / 12d);

            // KCIBUL 09-20-2010:
            double cu = u-m;
            return (cu - (0.5 * (cu<0?-1:1))) / sigma;
        }
        else {
            Histogram<Double> distribution = getDistribution(n1, n2);
            return (u - distribution.getMean()) / distribution.getStandardDeviation();

        }
    }

    /** Finds or calculates the median value of a sorted array of double. */
    public double median(final int[] data) {
        final int len = data.length;
        final int mid = len / 2;
        if (data.length % 2 == 0) {
            return (data[mid] + data[mid-1]) / 2d;
        }
        else {
            return data[mid];
        }
    }

    /** Constructs a new rank sum test with the given data. */
    public Result test(final int[] series1, final int[] series2) {
        final int n1 = series1.length;
        final int n2 = series2.length;

        final double u = calculateU(series1,  series2);
        final double z = calculateZ(u, n1, n2);
        double p;

        if (n1 >= this.minimumNormalN && n2 >= this.minimumNormalN) {
            // KCIBUL: 09-20-2010:
            p = 2 * NORMAL.cdf(NORMAL_MEAN - Math.abs(z) * NORMAL_SD);
        }
        else {
            Histogram<Double> distribution = getDistribution(n1, n2);
            double tmp = getCumulativeProbability(distribution, u);

            // KCIBUL: 09-17-2010: fixed bug where if tmp > 0.5, it was being set to tmp - 0.5!!
            if (tmp > 0.5) tmp = 1-tmp;

            // KCIBUL: 09-20-2010: don't need to take 2x since we chose the smallest u in the distribution
            p = tmp;
        }

        return new Result(u, z, p, Math.abs(median(series1) - median(series2)));
    }

    public double getCumulativeProbability(final Histogram<Double> histo, final double v) {
        double count = 0;
        double total = 0;

        for (final Histogram.Bin bin : histo.values()) {
            final double binValue = bin.getIdValue();
            if (binValue <= v) count += bin.getValue();
            total += bin.getValue();
        }

        return count / total;
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // Code to generate and cache distributions of U for pairs of data series with small
    // numbers of elements
    /////////////////////////////////////////////////////////////////////////////////////

    /** Key class to hold the size of the two series. */
    private static class Key {
        final int n1, n2;

        private Key(int n1, int n2) {
            this.n1 = Math.min(n1, n2);
            this.n2 = Math.max(n1, n2);
        }

        @Override public boolean equals(Object o) {
            if (o == null || getClass() != o.getClass()) return false;

            Key that = (Key) o;
            return (this.n1 == that.n1 && this.n2 == that.n2);
        }

        @Override public int hashCode() {
            int result = n1;
            return 31 * result + n2;
        }
    }


    Histogram<Double> getDistribution(int n1, int n2) {
        Key key = new Key(n1,n2);
        Histogram<Double> distribution = SMALL_DISTRIBUTIONS.get(key);
        if (distribution == null) {
            distribution = generateDistribution(n1 , n2);
            SMALL_DISTRIBUTIONS.put(key, distribution);
        }

        return distribution;
    }

    private Histogram<Double> generateDistribution(int n1, int n2) {
        final Histogram<Double> histo = new Histogram<Double>();
        final int[] series1 = new int[n1];
        final int[] series2 = new int[n2];
        final AbstractDistribution normal = this.modelDistribution;

        for (int i=0; i< valuesInCustomDistribution; ++i) {
            for (int j=0; j<n1; ++j) series1[j] = normal.nextInt();
            for (int j=0; j<n2; ++j) series2[j] = normal.nextInt();
            histo.increment(calculateU(series1, series2));
        }

        return histo;
    }

    public static void main(String[] args) {
        double TOL = 0.01d;
        RankSumTest rst = new RankSumTest();
        int[] a;
        int[] r;

        // generate 5 test cases
//        Random rand = new Random();
//        int minSize = 20;
//        int maxSize = 100;
//        int maxInt = 100;
//        for(int i=0;i<5;i++) {
//
//            int len = minSize + rand.nextInt(maxSize-minSize);
//            List<Integer> l = new ArrayList<Integer>();
//            for(int j=0; j<len; j++) {
//                l.add(rand.nextInt(maxInt));
//            }
//            System.out.println("a = [" + StringUtils.join(l,",") + "];");
//
//            len = minSize + rand.nextInt(maxSize-minSize);
//            l.clear();
//            for(int j=0; j<len; j++) {
//                l.add(rand.nextInt(maxInt));
//            }
//            System.out.println("r = [" + StringUtils.join(l,",") + "];");
//
//        }

        double matlabResult;
        double p;

        // 5 small
        a = new int[]{35,83,55,80,52,40,46};
        r = new int[]{74,87,89,0,70,55};
        matlabResult = 0.3112;
        p = rst.test(a,r).getP();
        if (Math.abs(p - matlabResult) > TOL) { throw new RuntimeException("not equal"); }

        a = new int[]{41,75,70,62,26,39,54,25,25,51,35,11,28,19,62,56};
        r = new int[]{8,41,75,61,93,70,61,84,4};
        matlabResult = 0.2124;
        p = rst.test(a,r).getP();
        if (Math.abs(p - matlabResult) > TOL) { throw new RuntimeException("not equal"); }

        a = new int[]{38,47,29,36,55,2,21,90,30,97,61,29,77,4,63,42};
        r = new int[]{71,74,89,59,38};
        matlabResult = 0.1369;
        p = rst.test(a,r).getP();
        if (Math.abs(p - matlabResult) > TOL) { throw new RuntimeException("not equal"); }

//        a = new int[]{73,80,2,87,54,43,28,55,35,45};
//        r = new int[]{34,47,17,12,86,24,46};
//        matlabResult = 0.3148;
//        p = rst.test(a,r).getP();
//        Assert.assertEquals("not equal", matlabResult, p, TOL);

        a = new int[]{65,55,37,15,94,46,31,96,97,97,92,38};
        r = new int[]{21,76,41,66,22,58,90};
        matlabResult = 0.4714;
        p = rst.test(a,r).getP();
        if (Math.abs(p - matlabResult) > TOL) { throw new RuntimeException("not equal"); }





        // 5 longer
        a = new int[]{9,21,16,52,91,14,1,44,58,17,28,70,82,9,80,97,36,80,96,20,85,40,59,1,6,18,81,30,14,3,39,3,52,2,28,55,66,78,15,50,83,85,97,52,53};
        r = new int[]{81,93,60,58,63,58,80,46,18,65,53,49,2,37,48,84,76,95,10,13,77,77,97,67,83,25,1,31,92,9,12,21,12,35,17,23,86,4,88,74,95,72,7};
        matlabResult = 0.3784;
        p = rst.test(a,r).getP();
        if (Math.abs(p - matlabResult) > TOL) { throw new RuntimeException("not equal"); }

        a = new int[]{67,29,33,29,94,75,67,34,26,0,68,91,3,16,75,65,44,13,92,69,18,82,60,84,67,34,14,9,55,92,74,96,43,86,70,85,48,98,16,77,4,29,78,40,92,15,65,69,57,56,9,46,17,19,48,37,78,4,86,81,83,43,28,53,91,39,68,88,94,53,64,9,42,12,33,49,43,27,24,47,0,14,47,40,42,17,49,38};
        r = new int[]{33,75,95,45,1,20,30,78,70,62,50,75,71,75,67,69,77,71,62,70,35,58,37,65,26,86,32,77,32};
        matlabResult = 0.2352;
        p = rst.test(a,r).getP();
        if (Math.abs(p - matlabResult) > TOL) { throw new RuntimeException("not equal"); }

        a = new int[]{26,69,31,51,92,71,81,90,37,77,30,67,11,33,29,98,70,65,24,76,68,47,8,73,9,87,78,65,4,19,89,39,18,96,20};
        r = new int[]{55,47,75,95,23,19,50,46,72,55,2,36,15,96,81,46,33,60,13,5,31,27,24,53,94,94,31,96,45,10,9,47,24,95,0,58,95,45,15,48,0,49,90,33,82,97,48,7,95,63,70,96,65,52,86,2,80,52,11,34,40,38,7,88,91,93,11,23,27,59,70,25,23,75,2,17,62,62,20,72,72,94,12,16,16,46,54,0,8,55,52,26,51,21,76,78};
        matlabResult = 0.3681;
        p = rst.test(a,r).getP();
        if (Math.abs(p - matlabResult) > TOL) { throw new RuntimeException("not equal"); }

        a = new int[]{19,16,51,4,16,93,32,89,81,22,9,52,71,63,33,78,55,34,60,3,27,69,14,21,56,94,55,89,86,41,6,36,95,12,77,1,28,57,71,70,51};
        r = new int[]{28,20,65,83,78,36,84,93,9,40,4,95,85,8,64,66,25,43,62,30,68,66,75,18,73,55,8,62};
        matlabResult = 0.5251;
        p = rst.test(a,r).getP();
        if (Math.abs(p - matlabResult) > TOL) { throw new RuntimeException("not equal"); }

        a = new int[]{62,77,38,22,81,64,62,30,19,11,60,5,56,84,24,6,15,53,35,16,1,36,26,28,81,34,31,3,66,12,92,32,9,63,60,0,0,55,58,79,85,87,11,42,35,49,32,98};
        r = new int[]{49,42,33,22,57,99,38,83,35,34,48,37,27,95,82,52,61,16,69,5,83,56,39,23,51,31,84,68};
        matlabResult = 0.1765;
        p = rst.test(a,r).getP();
        if (Math.abs(p - matlabResult) > TOL) { throw new RuntimeException("not equal"); }






        // specific test cases we've had trouble with
        a = new int[]{2,5,8,19,21,21,31,40,43,45,71,72};
        r = new int[]{0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,6,6,6,7,7,7,7,8,8,8,8,8,8,8,8,8,9,10,10,10,10,10,10,11,11,11,12,12,12,13,14,14,14,14,15,15,15,15,16,17,17,17,17,18,18,19,19,19,19,20,20,20,21,21,21,21,21,21,22,22,22,22,23,23,23,23,24,24,24,24,24,24,24,25,25,25,25,25,25,25,25,25,26,26,26,26,27,27,27,27,27,27,27,28,28,28,29,29,29,29,29,29,29,30,31,31,31,31,31,33,34,35,36,36,37,39,39,40,40,42,43,44,44,44,45,45,45,46,46,46,47,47,48,48,49,49,49,49,50,51,51,51,52,52,52,53,53,55,57,58,60,61,61,62,62,62,63,64,65,66,67,68,71,71,71,72,75};
        matlabResult = 0.4969;
        p = rst.test(a,r).getP();
        if (Math.abs(p - matlabResult) > TOL) { throw new RuntimeException("not equal"); }

        
    }
}
