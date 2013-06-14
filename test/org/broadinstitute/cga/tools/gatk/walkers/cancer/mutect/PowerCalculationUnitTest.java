package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import org.apache.commons.math.MathException;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Unit tests for MuTect
 */
public class PowerCalculationUnitTest {
    @Test
    public void testTumor() throws MathException {
        double Q30_EPS = Math.pow(10, (-30/10)) / 3d;
        double LOD = 6.3d;
        double wiggle = 0.00001;

        // test log likelihood calculation
        Assert.assertEquals(TumorPowerCalculator.calculateLogLikelihood(50, 5, 0.001, 0.1), -7.059164, wiggle, "Failed!");
        Assert.assertEquals(TumorPowerCalculator.calculateLogLikelihood(40, 5, 0.001, 0.1), -6.597727, wiggle, "Failed!");
        Assert.assertEquals(TumorPowerCalculator.calculateLogLikelihood(40, 10, 0.001, 0.1), -11.34971, wiggle, "Failed!");
        Assert.assertEquals(TumorPowerCalculator.calculateLogLikelihood(40, 10, 0.01, 0.1), -11.15482, wiggle, "Failed!");
        Assert.assertEquals(TumorPowerCalculator.calculateLogLikelihood(40, 10, 0.01, 0.2), -9.866713, wiggle, "Failed!");


        Assert.assertEquals(testTumorCalculatePower(50, Q30_EPS, LOD, 0.15, 0), 0.975342, wiggle, "Failed!");

        Assert.assertEquals(testTumorCalculatePower(20, Q30_EPS, LOD, 0.15, 0), 0.6364833, wiggle, "Failed!");
        Assert.assertEquals(testTumorCalculatePower(10, Q30_EPS, LOD, 0.15, 0), 0.3167154, wiggle, "Failed!");
        Assert.assertEquals(testTumorCalculatePower(50, Q30_EPS, LOD, 0.35, 0), 0.9999994, wiggle, "Failed!");
        Assert.assertEquals(testTumorCalculatePower(50, Q30_EPS, LOD, 0.05, 0), 0.3893047, wiggle, "Failed!");
        Assert.assertEquals(testTumorCalculatePower(50, Q30_EPS, 4.3, 0.05, 0), 0.6068266, wiggle, "Failed!");
        Assert.assertEquals(testTumorCalculatePower(50, Q30_EPS, 4.3, 0.05, 0.02), 0.0039135, wiggle, "Failed!");

        // make an instance to test caching
        TumorPowerCalculator pc = new TumorPowerCalculator(Q30_EPS, LOD, 0);
        Assert.assertEquals(pc.cachingPowerCalculation(50, 0.15), 0.975342, wiggle, "Failed!");
        Assert.assertEquals(pc.cachingPowerCalculation(50, 0.15), 0.975342, wiggle, "Failed!");
        Assert.assertEquals(pc.cachingPowerCalculation(50, 0.15), 0.975342, wiggle, "Failed!");

    }

    private static double testTumorCalculatePower(int n, double eps, double fdr, double delta, double contam) throws MathException {
        double power = TumorPowerCalculator.calculatePower(n, eps, fdr, delta, contam, true);
        System.out.println("Depth: " + n + " EPS: " + eps + " FDR: " + fdr + " Delta: " + delta + " --> Power: " + power);
        return power;
    }



    @Test
    public void testNormal() throws MathException {
        double Q30_EPS = Math.pow(10, (-30/10)) / 3d;
        double LOD = 2.3d;
        double wiggle = 0.000001;

        Assert.assertEquals(testNormalCalculatePower(50, Q30_EPS, LOD), 1.0000000, wiggle, "Failed!");
        Assert.assertEquals(testNormalCalculatePower(10, Q30_EPS, LOD), 0.9973492, wiggle, "Failed!");
        Assert.assertEquals(testNormalCalculatePower( 8, Q30_EPS, LOD), 0.9974184, wiggle, "Failed!");
        Assert.assertEquals(testNormalCalculatePower( 5, Q30_EPS, LOD), 0.0000000, wiggle, "Failed!");
        Assert.assertEquals(testNormalCalculatePower(10, Math.pow(10, (-20/10)) / 3d, LOD), 0.9762534, wiggle, "Failed!");


        // make an instance to test caching
        NormalPowerCalculator pc = new NormalPowerCalculator(Q30_EPS, LOD, true);
        Assert.assertEquals(pc.cachingPowerCalculation(10), 0.9973492, wiggle, "Failed!");
        Assert.assertEquals(pc.cachingPowerCalculation(10), 0.9973492, wiggle, "Failed!");
        Assert.assertEquals(pc.cachingPowerCalculation(10), 0.9973492, wiggle, "Failed!");

    }

    private static double testNormalCalculatePower(int n, double eps, double threshold) throws MathException {
        double power = NormalPowerCalculator.calculatePower(n, eps, threshold, true);
        System.out.println("Depth: " + n + " EPS: " + eps + " Threshold: " + threshold + " --> Power: " + power);
        return power;
    }



}

