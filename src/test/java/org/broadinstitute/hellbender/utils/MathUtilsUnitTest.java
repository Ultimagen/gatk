package org.broadinstitute.hellbender.utils;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.internal.junit.ArrayAsserts;

import java.util.*;
import java.util.function.DoublePredicate;
import java.util.function.DoubleUnaryOperator;
import java.util.function.IntPredicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static java.lang.Math.exp;
import static java.lang.Math.log10;
import static org.broadinstitute.hellbender.utils.MathUtils.*;

/**
 * Basic unit test for MathUtils
 */
public final class MathUtilsUnitTest extends GATKBaseTest {

    private static final Logger logger = LogManager.getLogger(MathUtilsUnitTest.class);

    /**
     * Tests that we correctly compute mean and standard deviation from a stream of numbers
     */
    @Test
    public void testRunningAverage() {
        logger.warn("Executing testRunningAverage");

        int[] numbers = {1, 2, 4, 5, 3, 128, 25678, -24};
        MathUtils.RunningAverage r = new MathUtils.RunningAverage();

        for (final double b : numbers)
            r.add(b);

        Assert.assertEquals((long) numbers.length, r.observationCount());
        Assert.assertTrue(r.mean() - 3224.625 < 2e-10);
        Assert.assertTrue(r.stddev() - 9072.6515881128 < 2e-10);
    }

    @Test(dataProvider = "log10OneMinusPow10Data")
    public void testLog10OneMinusPow10(final double x, final double expected) {
        final double actual = MathUtils.log10OneMinusPow10(x);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual,expected,1E-9);
    }

    @Test(dataProvider = "log1mexpData")
    public void testLog1mexp(final double x, final double expected) {
        final double actual = NaturalLogUtils.log1mexp(x);
        if (Double.isNaN(expected))
            Assert.assertTrue(Double.isNaN(actual));
        else
            Assert.assertEquals(actual,expected,1E-9);
    }

    @DataProvider(name = "log10OneMinusPow10Data")
    public Iterator<Object[]> log10OneMinusPow10Data() {

        final double[] inValues = new double[] { Double.NaN, 10, 1, 0, -1, -3, -10, -30, -100, -300, -1000, -3000 };
        return new Iterator<Object[]>() {

            private int i = 0;

            @Override
            public boolean hasNext() {
                return i < inValues.length;

            }

            @Override
            public Object[] next() {
                final double input = inValues[i++];
                final double output = log10(1 - Math.pow(10, input));
                return new Object[] { input, output };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @DataProvider(name = "log1mexpData")
    public Iterator<Object[]> log1mexpData() {

        final double[] inValues = new double[] { Double.NaN, 10, 1, 0, -1, -3, -10, -30, -100, -300, -1000, -3000 };
        return new Iterator<Object[]>() {

            private int i = 0;

            @Override
            public boolean hasNext() {
                return i < inValues.length;

            }

            @Override
            public Object[] next() {
                final double input = inValues[i++];
                final double output = Math.log(1 - exp(input));
                return new Object[] { input, output };
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @Test
    public void testGoodProbability(){
        for (final double good : Arrays.asList(0.1, 1.0)) {
            Assert.assertTrue(MathUtils.isValidProbability(good));
            Assert.assertTrue(MathUtils.isValidLog10Probability(log10(good)));
        }

        Assert.assertTrue(MathUtils.isValidProbability(0.0));
        Assert.assertTrue(MathUtils.isValidLog10Probability(Double.NEGATIVE_INFINITY));

        for (final double bad : Arrays.asList(-1.0, 2.0, Double.NaN, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY)) {
            Assert.assertFalse(MathUtils.isValidProbability(bad));
        }
    }

    @Test
    public void testLog10OneMinusX(){
        Assert.assertEquals(-0.04575749056, MathUtils.log10OneMinusX(0.1), 1e-6); //result from wolphram alpha
        Assert.assertEquals(0.0, MathUtils.log10OneMinusX(0.0), 1e-6);
        Assert.assertEquals(Double.NEGATIVE_INFINITY, MathUtils.log10OneMinusX(1.0), 1e-6);
    }

    @Test
    public void testBinomialProbability() {
        // results from Wolfram Alpha
        Assert.assertEquals(MathUtils.binomialProbability(3, 2, 0.5), 0.375, 0.0001);
        Assert.assertEquals(MathUtils.binomialProbability(100, 10, 0.5), 1.365543e-17, 1e-18);
        Assert.assertEquals(MathUtils.binomialProbability(217, 73, 0.02), 4.521904e-67, 1e-68);
        Assert.assertEquals(MathUtils.binomialProbability(300, 100, 0.02), 9.27097e-91, 1e-92);
        Assert.assertEquals(MathUtils.binomialProbability(300, 150, 0.98), 6.462892e-168, 1e-169);
        Assert.assertEquals(MathUtils.binomialProbability(300, 120, 0.98), 3.090054e-221, 1e-222);
        Assert.assertEquals(MathUtils.binomialProbability(300, 112, 0.98), 2.34763e-236, 1e-237);
    }

    @Test
    public void testApproximateLogSumLog() {
        final double requiredPrecision = 1E-4;

        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{0.0, 0.0, 0.0}), log10(3), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{0.0, 0.0, 0.0}, 0), 0.0, requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{0.0, 0.0, 0.0}, 3), log10(3), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{0.0, 0.0, 0.0}, 2), log10(2), requiredPrecision);
        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{0.0, 0.0, 0.0}, 1), 0.0, requiredPrecision);

        Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY}), Double.NEGATIVE_INFINITY);

        final Random random = new Random(13);
        for (int j = 0; j < 5; j++) {
            for (int i = 0; i < 5; i++) {
                final double a = (1 + 3 * j) * random.nextGaussian();
                final double b = (1 + 3 * j) * random.nextGaussian();
                final double c = (1 + 3 * j) * random.nextGaussian();

                Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{a}), a, requiredPrecision);
                Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{a, Double.NEGATIVE_INFINITY}), a, requiredPrecision);
                Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{a, b}), log10(exp10(a) + exp10(b)), requiredPrecision);
                Assert.assertEquals(MathUtils.approximateLog10SumLog10(a, b), log10(exp10(a) + exp10(b)), requiredPrecision);
                Assert.assertEquals(MathUtils.approximateLog10SumLog10(b, a), log10(exp10(a) + exp10(b)), requiredPrecision);
                Assert.assertEquals(MathUtils.approximateLog10SumLog10(new double[]{a, b, c}), log10(exp10(a) + exp10(b) + exp10(c)), requiredPrecision);
                Assert.assertEquals(MathUtils.approximateLog10SumLog10(a, b, c), log10(exp10(a) + exp10(b) + exp10(c)), requiredPrecision);
            }
        }
    }

    @Test
    public void testNormalDistribution() {
        final double requiredPrecision = 1E-10;

        for( final double mu : new double[]{-5.0, -3.2, -1.5, 0.0, 1.2, 3.0, 5.8977} ) {
            for( final double sigma : new double[]{1.2, 3.0, 5.8977} ) {
                for( final double x : new double[]{-5.0, -3.2, -1.5, 0.0, 1.2, 3.0, 5.8977} ) {
                    // TODO: GATK3 uses the cern.jet.random Normal implementation to verify this;
                    // since it was only used for test verification, rather than introducing a new
                    // GATK4 dependency we just use the apache NormalDistribution
                    final NormalDistribution n = new NormalDistribution(mu, sigma);
                    Assert.assertEquals(n.density(x), MathUtils.normalDistribution(mu, sigma, x), requiredPrecision);
                    Assert.assertEquals(Math.log10(n.density(x)), MathUtils.normalDistributionLog10(mu, sigma, x), requiredPrecision);
                }
            }
        }
    }

    private static double exp10(final double x){
        return Math.pow(10.0, x);
    }

    @Test
    public void testLogSumLog() {
        final double requiredPrecision = 1E-14;

        Assert.assertEquals(MathUtils.log10SumLog10(new double[]{0.0, 0.0, 0.0}), log10(3), requiredPrecision);
        Assert.assertEquals(MathUtils.log10SumLog10(new double[]{0.0, 0.0, 0.0}, 0), log10(3), requiredPrecision);
        Assert.assertEquals(MathUtils.log10SumLog10(new double[]{0.0, 0.0, 0.0}, 0, 3), log10(3), requiredPrecision);
        Assert.assertEquals(MathUtils.log10SumLog10(new double[]{0.0, 0.0, 0.0}, 0, 2), log10(2), requiredPrecision);
        Assert.assertEquals(MathUtils.log10SumLog10(new double[]{0.0, 0.0, 0.0}, 0, 1), 0.0, requiredPrecision);

        final Random random = new Random(13);
        for (int j = 0; j < 5; j++) {
            for (int i = 0; i < 5; i++) {
                final double a = (1 + 3 * j) * random.nextGaussian();
                final double b = (1 + 3 * j) * random.nextGaussian();
                final double c = (1 + 3 * j) * random.nextGaussian();

                Assert.assertEquals(MathUtils.log10SumLog10(new double[]{a}), a, requiredPrecision);
                Assert.assertEquals(MathUtils.log10SumLog10(new double[]{a, Double.NEGATIVE_INFINITY}), a, requiredPrecision);
                Assert.assertEquals(MathUtils.log10SumLog10(new double[]{a, b}), log10(exp10(a) + exp10(b)), requiredPrecision);
                Assert.assertEquals(MathUtils.log10SumLog10(new double[]{a, b, c}), log10(exp10(a) + exp10(b) + exp10(c)), requiredPrecision);
            }
        }
    }

    @Test
    public void testlogsumLogEdgeCases(){
        final double requiredPrecision = 1E-14;

        Assert.assertEquals(MathUtils.log10SumLog10(new double[]{3.0, 2.0, 1.0}, 2, 1), Double.NEGATIVE_INFINITY);
        Assert.assertEquals(MathUtils.log10SumLog10(new double[]{Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY}), Double.NEGATIVE_INFINITY);
        Assert.assertEquals(MathUtils.log10SumLog10(new double[]{3.0, 2.0, 1.0}), log10(exp10(3.0) + exp10(2.0) + exp10(1.0)), requiredPrecision);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testlogsumLogErrorNaN(){
        MathUtils.log10SumLog10(new double[]{Double.NaN, 1.0});
    }

    @Test
    public void testNormalize(){
        final double error = 1e-6;

        final double[] normalized      = MathUtils.normalizeFromLog10ToLinearSpace(new double[]{log10(3.0), log10(2.0), log10(1.0)});
        final double[] normalizedLog10 = MathUtils.normalizeLog10(new double[]{log10(3.0), log10(2.0), log10(1.0)});
        final double[] normalizedLogInLog10 = MathUtils.scaleLogSpaceArrayForNumericalStability(new double[]{3.0, 2.0, 1.0});
        final double[] normalizedExpected    = {3.0/6.0, 2.0/6.0, 1.0/6.0};  //sum of those is 1
        final double[] normalizedLog10Expected = {log10(3.0 / 6.0), log10(2.0 / 6.0), log10(1.0 / 6.0)}; //sum of 10^ of those is 1
        final double[] normalizedLogInLogExpected = {0.0, -1.0, -2.0};

        assertEqualsDoubleArray(normalizedExpected, normalized, error);
        assertEqualsDoubleArray(normalizedLog10Expected, normalizedLog10, error);
        assertEqualsDoubleArray(normalizedLogInLogExpected, normalizedLogInLog10, error);
    }

    @Test
    public void testSum() {
        double[] doubleTest = {-1,0,1,2,3};
        long[] longTest = {-1,0,1,2,3};
        int[] intTest = {-1,0,1,2,3};
        Assert.assertEquals(MathUtils.sum(doubleTest), 5.0);
        Assert.assertEquals(MathUtils.sum(longTest), 5);
        Assert.assertEquals(MathUtils.sum(intTest), 5);
    }

    @Test
    public void testSumLog10() {
        double[] xLog10 = {log10(1.0/6), log10(2.0/6), log10(3.0/6)};
        Assert.assertEquals(MathUtils.sumLog10(xLog10), 1.0);
    }

    @Test
    public void testSumRange() {
        double[] doubleTest = {-1,0,1,2,3};
        Assert.assertEquals(MathUtils.sum(doubleTest, 1, 4), 3.0);
    }

    @Test
    public void testSumCollection() {
        final List<Integer> list = Arrays.asList(1,2,3,4,5);
        Assert.assertEquals(MathUtils.sumDoubleFunction(list, n -> 0.5 * n), 7.5);
    }

    @Test
    public void testCompareDoubles(){
        Assert.assertEquals(MathUtils.compareDoubles(0.1, 0.2), 1);
        Assert.assertEquals(MathUtils.compareDoubles(0.2, 0.1), -1);
        Assert.assertEquals(MathUtils.compareDoubles(0.1, 0.1), 0);
        Assert.assertEquals(MathUtils.compareDoubles(0.1, 0.1 + (1e-7)), 0);
        Assert.assertEquals(MathUtils.compareDoubles(0.1, 0.1 + (1e-7), 1e-8), 1);
    }

    @Test
    public void testNormalizeFromReal(){
        final double error = 1e-6;
        final double[] actual = MathUtils.normalizeSumToOne(new double[]{1.0, 2.0, 3.0});
        final double[] expected =  {1.0/6.0, 2.0/6.0, 3.0/6.0};
        for (int i = 0; i < actual.length; i++){
            Assert.assertEquals(expected[i], actual[i], error);
        }
    }

    @Test
    public void testLogLog10conversions(){
        final double error = 1e-6;
        for (final double x : new double[]{Math.E, 10.0, 0.5, 1, 1.5, 100.0}) {
            final double logX = Math.log(x);
            final double log10X = Math.log10(x);
            Assert.assertEquals(log10ToLog(log10X), logX, error, "log10ToLog");
            Assert.assertEquals(logToLog10(logX), log10X, error, "logToLog10");
            Assert.assertEquals(logToLog10(log10ToLog(log10X)), log10X, error, "log10->log->log10");
            Assert.assertEquals(log10ToLog(logToLog10(logX)), logX, error, "log->log10->log");
        }
    }


    @DataProvider(name = "rounding")
    public Object[][] makeRounding() {
        List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{3.1415926, 8, 3.1415926});
        tests.add(new Object[]{3.1415926, 7, 3.1415926});
        tests.add(new Object[]{3.1415926, 6, 3.141593});
        tests.add(new Object[]{3.1415926, 5, 3.14159});
        tests.add(new Object[]{3.1415926, 4, 3.1416});
        tests.add(new Object[]{3.1415926, 3, 3.142});
        tests.add(new Object[]{3.1415926, 2, 3.14});
        tests.add(new Object[]{3.1415926, 1, 3.1});

        tests.add(new Object[]{1.025, 2, 1.03});
        tests.add(new Object[]{34.356845398701736,  4, 34.3568});
        tests.add(new Object[]{0.43452380952380953, 2, 0.43});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "rounding")
    public void testRounding(final double in, final int n, final double out) {
        Assert.assertEquals(roundToNDecimalPlaces(in, n), out);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testRoundingError() {
        roundToNDecimalPlaces(1.1234, -1);
    }

    /**
     * This is not really a test but a demo that adding doubles can have different results depending on the order.
     * You can remedy that my using multipliers but that too can lead to different results depending on the multiplier.
     */
    @Test
    public void testAddDoubles() throws Exception {
        final double[] ds =
                {0.125,
                 0.125,
                 0.125, 0.125, 0.125, 0.125,
                        0.125, 0.14285714285714285, 0.14285714285714285, 0.14285714285714285,
                        0.14285714285714285, 0.14285714285714285, 0.14285714285714285, 0.14285714285714285,
                        0.125, 0.125, 0.125, 0.125,
                        0.125, 0.125, 0.125, 0.125};

        final double[] ds_reordered = {
        0.125,
        0.125,
        0.125,
        0.125,
        0.125,
        0.125,
        0.125,
        0.125,
        0.125,
        0.125,
        0.125,
        0.125,
        0.125,
        0.125,
        0.125,
        0.14285714285714285,
        0.14285714285714285,
        0.14285714285714285,
        0.14285714285714285,
        0.14285714285714285,
        0.14285714285714285,
        0.14285714285714285};

        //Use different multipliers, see that the results are all different
        double sumNoMult = 0.0;
        for (int i = 0; i < ds.length; i++) {
            sumNoMult += ds[i];
        }
        double sumNoMult_reordered = 0.0;
        for (int i = 0; i < ds_reordered.length; i++) {
            sumNoMult_reordered += ds[i];
        }

        double sumMult1000 = 0.0;
        for (int i = 0; i < ds.length; i++) {
            sumMult1000 += ds[i] * 1000.0;
        }
        sumMult1000 /= 1000.0;

        double sumMult10000 = 0.0;
        for (int i = 0; i < ds.length; i++) {
            sumMult10000 += ds[i] * 10000.0;
        }
        sumMult10000 /= 10000.0;

        double sumMult100000 = 0.0;
        for (int i = 0; i < ds.length; i++) {
            sumMult100000 += ds[i] * 100000.0;
        }
        sumMult100000 /= 100000.0;

        double sumMult100000_reordered = 0.0;
        for (int i = 0; i < ds_reordered.length; i++) {
            sumMult100000_reordered += ds_reordered[i] * 100000.0;
        }
        sumMult100000_reordered /= 100000.0;

        double sumMult1000_reordered = 0.0;
        for (int i = 0; i < ds_reordered.length; i++) {
            sumMult1000_reordered += ds_reordered[i] * 1000.0;
        }
        sumMult1000_reordered /= 1000.0;

        //They are all different
        Assert.assertNotEquals(sumNoMult, sumMult1000);
        Assert.assertNotEquals(sumNoMult, sumMult10000);
        Assert.assertNotEquals(sumMult1000, sumMult10000);

        Assert.assertNotEquals(sumMult1000, sumMult1000_reordered); //not equal to itself when reordered

        //But these ones is the same ordered or not (though not equal to each other)
        Assert.assertEquals(sumNoMult, sumNoMult_reordered, sumNoMult + " vs " + sumNoMult_reordered);
        Assert.assertEquals(sumMult100000_reordered, sumMult100000);
        Assert.assertNotEquals(sumNoMult, sumMult100000);
    }


    @DataProvider(name="numbersForMedian")
    public Object[][] getNumbersForMedian(){
        return new Object[][] {
                {Arrays.asList(0), 0},
                {Arrays.asList(0,1, 2), 1},
                {Arrays.asList(1, 0, 2), 1},
                {Arrays.asList(0, 1), .5},
                {Arrays.asList(-10.0, 0, 10.5), 0},
                {Arrays.asList(1, 1, 5, 2, 1, 5, 1, 9),1.5},
        };
    }

    @Test(dataProvider = "numbersForMedian")
    public <T extends Number & Comparable<T>> void testMedian(List<T> values, double expected){
        Assert.assertEquals(MathUtils.median(values), expected);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMedianOfEmptyList(){
        Collection<Integer> empty = Collections.emptyList();
        MathUtils.median(empty);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMedianOfNullList(){
        MathUtils.median((List<Integer>) null);
    }

    @Test
    public void testLog10sumLog10() {
        final double requiredPrecision = 1E-14;

        final double log3 = 0.477121254719662;
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{0.0, 0.0, 0.0}), log3, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {0.0, 0.0, 0.0}, 0), log3, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[]{0.0, 0.0, 0.0}, 0, 3), log3, requiredPrecision);

        final double log2 = 0.301029995663981;
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {0.0, 0.0, 0.0}, 0, 2), log2, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {0.0, 0.0, 0.0}, 0, 1), 0.0, requiredPrecision);

        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {0.0}), 0.0, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-5.15}), -5.15, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {130.0}), 130.0, requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-0.145}), -0.145, requiredPrecision);

        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {0.0, 0.0}), Math.log10(Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-1.0, 0.0}), Math.log10(Math.pow(10.0, -1.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {0.0, -1.0}), Math.log10(Math.pow(10.0, 0.0) + Math.pow(10.0, -1.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-2.2, -3.5}), Math.log10(Math.pow(10.0, -2.2) + Math.pow(10.0, -3.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-1.0, -7.1}), Math.log10(Math.pow(10.0, -1.0) + Math.pow(10.0, -7.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {5.0, 6.2}), Math.log10(Math.pow(10.0, 5.0) + Math.pow(10.0, 6.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {38.1, 16.2}), Math.log10(Math.pow(10.0, 38.1) + Math.pow(10.0, 16.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-38.1, 6.2}), Math.log10(Math.pow(10.0, -38.1) + Math.pow(10.0, 6.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-19.1, -37.1}), Math.log10(Math.pow(10.0, -19.1) + Math.pow(10.0, -37.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-29.1, -27.6}), Math.log10(Math.pow(10.0, -29.1) + Math.pow(10.0, -27.6)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-0.12345, -0.23456}), Math.log10(Math.pow(10.0, -0.12345) + Math.pow(10.0, -0.23456)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-15.7654, -17.0101}), Math.log10(Math.pow(10.0, -15.7654) + Math.pow(10.0, -17.0101)), requiredPrecision);

        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {0.0, 0.0, 0.0}), Math.log10(Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-1.0, 0.0, 0.0}), Math.log10(Math.pow(10.0, -1.0) + Math.pow(10.0, 0.0) + Math.pow(10.0, 0.0)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {0.0, -1.0, -2.5}), Math.log10(Math.pow(10.0, 0.0) + Math.pow(10.0, -1.0) + Math.pow(10.0, -2.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-2.2, -3.5, -1.1}), Math.log10(Math.pow(10.0, -2.2) + Math.pow(10.0, -3.5) + Math.pow(10.0, -1.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-1.0, -7.1, 0.5}), Math.log10(Math.pow(10.0, -1.0) + Math.pow(10.0, -7.1) + Math.pow(10.0, 0.5)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {5.0, 6.2, 1.3}), Math.log10(Math.pow(10.0, 5.0) + Math.pow(10.0, 6.2) + Math.pow(10.0, 1.3)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {38.1, 16.2, 18.1}), Math.log10(Math.pow(10.0, 38.1) + Math.pow(10.0, 16.2) + Math.pow(10.0, 18.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-38.1, 6.2, 26.6}), Math.log10(Math.pow(10.0, -38.1) + Math.pow(10.0, 6.2) + Math.pow(10.0, 26.6)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-19.1, -37.1, -45.1}), Math.log10(Math.pow(10.0, -19.1) + Math.pow(10.0, -37.1) + Math.pow(10.0, -45.1)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-29.1, -27.6, -26.2}), Math.log10(Math.pow(10.0, -29.1) + Math.pow(10.0, -27.6) + Math.pow(10.0, -26.2)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-0.12345, -0.23456, -0.34567}), Math.log10(Math.pow(10.0, -0.12345) + Math.pow(10.0, -0.23456) + Math.pow(10.0, -0.34567)), requiredPrecision);
        Assert.assertEquals(MathUtils.log10sumLog10(new double[] {-15.7654, -17.0101, -17.9341}), Math.log10(Math.pow(10.0, -15.7654) + Math.pow(10.0, -17.0101) + Math.pow(10.0, -17.9341)), requiredPrecision);

        // magnitude of the sum doesn't matter, so we can combinatorially test this via partitions of unity
        double[] mult_partitionFactor = new double[]{0.999,0.98,0.95,0.90,0.8,0.5,0.3,0.1,0.05,0.001};
        int[] n_partitions = new int[] {2,4,8,16,32,64,128,256,512,1028};
        for ( double alpha : mult_partitionFactor ) {
            double log_alpha = Math.log10(alpha);
            double log_oneMinusAlpha = Math.log10(1-alpha);
            for ( int npart : n_partitions ) {
                double[] multiplicative = new double[npart];
                double[] equal = new double[npart];
                double remaining_log = 0.0;  // realspace = 1
                for ( int i = 0 ; i < npart-1; i++ ) {
                    equal[i] = -Math.log10(npart);
                    double piece = remaining_log + log_alpha; // take a*remaining, leaving remaining-a*remaining = (1-a)*remaining
                    multiplicative[i] = piece;
                    remaining_log = remaining_log + log_oneMinusAlpha;
                }
                equal[npart-1] = -Math.log10(npart);
                multiplicative[npart-1] = remaining_log;
                Assert.assertEquals(MathUtils.log10sumLog10(equal),0.0,requiredPrecision);
                Assert.assertEquals(MathUtils.log10sumLog10(multiplicative),0.0,requiredPrecision,String.format("Did not sum to one: nPartitions=%d, alpha=%f",npart,alpha));
            }
        }
    }

    /**
     * The PartitionGenerator generates all of the partitions of a number n, e.g.
     * 5 + 0
     * 4 + 1
     * 3 + 2
     * 3 + 1 + 1
     * 2 + 2 + 1
     * 2 + 1 + 1 + 1
     * 1 + 1 + 1 + 1 + 1
     *
     * This is used to help enumerate the state space over which the Dirichlet-Multinomial is defined,
     * to ensure that the distribution function is properly implemented
     */
    class PartitionGenerator implements Iterator<List<Integer>> {
        // generate the partitions of an integer, each partition sorted numerically
        int n;
        List<Integer> a;

        int y;
        int k;
        int state;

        int x;
        int l;

        public PartitionGenerator(int n) {
            this.n = n;
            this.y = n - 1;
            this.k = 1;
            this.a = new ArrayList<>();
            for ( int i = 0; i < n; i++ ) {
                this.a.add(i);
            }
            this.state = 0;
        }

        @Override
        public void remove()  { /* do nothing */ }

        @Override
        public boolean hasNext() { return ! ( this.k == 0 && state == 0 ); }

        private String dataStr()  {
            return String.format("a = [%s]  k = %d  y = %d  state = %d  x = %d  l = %d",
                    Utils.join(",",a), k, y, state, x, l);
        }

        @Override
        public List<Integer> next() {
            if ( this.state == 0 ) {
                this.x = a.get(k-1)+1;
                k -= 1;
                this.state = 1;
            }

            if ( this.state == 1 ) {
                while ( 2 * x <= y ) {
                    this.a.set(k,x);
                    this.y -= x;
                    this.k++;
                }
                this.l = 1+this.k;
                this.state = 2;
            }

            if ( this.state == 2 ) {
                if ( x <= y ) {
                    this.a.set(k,x);
                    this.a.set(l,y);
                    x += 1;
                    y -= 1;
                    return this.a.subList(0, this.k + 2);
                } else {
                    this.state =3;
                }
            }

            if ( this.state == 3 ) {
                this.a.set(k,x+y);
                this.y = x + y - 1;
                this.state = 0;
                return a.subList(0, k + 1);
            }

            throw new IllegalStateException("Cannot get here");
        }

        public String toString() {
            final StringBuilder buf = new StringBuilder();
            buf.append("{ ");
            while ( hasNext() ) {
                buf.append("[");
                buf.append(Utils.join(",",next()));
                buf.append("],");
            }
            buf.deleteCharAt(buf.lastIndexOf(","));
            buf.append(" }");
            return buf.toString();
        }

    }

    /**
     * NextCounts is the enumerator over the state space of the multinomial dirichlet.
     *
     * It filters the partition of the total sum to only those with a number of terms
     * equal to the number of categories.
     *
     * It then generates all permutations of that partition.
     *
     * In so doing it enumerates over the full state space.
     */
    class NextCounts implements Iterator<int[]> {

        private PartitionGenerator partitioner;
        private int numCategories;
        private int[] next;

        public NextCounts(int numCategories, int totalCounts) {
            partitioner = new PartitionGenerator(totalCounts);
            this.numCategories = numCategories;
            next = nextFromPartitioner();
        }

        @Override
        public void remove() { /* do nothing */ }

        @Override
        public boolean hasNext() { return next != null; }

        @Override
        public int[] next() {
            int[] toReturn = clone(next);
            next = nextPermutation();
            if ( next == null ) {
                next = nextFromPartitioner();
            }

            return toReturn;
        }

        private int[] clone(int[] arr) {
            return Arrays.copyOf(arr, arr.length);
        }

        private int[] nextFromPartitioner() {
            if ( partitioner.hasNext() ) {
                List<Integer> nxt = partitioner.next();
                while ( partitioner.hasNext() && nxt.size() > numCategories ) {
                    nxt = partitioner.next();
                }

                if ( nxt.size() > numCategories ) {
                    return null;
                } else {
                    int[] buf = new int[numCategories];
                    for ( int idx = 0; idx < nxt.size(); idx++ ) {
                        buf[idx] = nxt.get(idx);
                    }
                    Arrays.sort(buf);
                    return buf;
                }
            }

            return null;
        }

        public int[] nextPermutation() {
            return MathUtilsUnitTest.nextPermutation(next);
        }

    }

    public static int[] nextPermutation(int[] next) {
        // the counts can swap among each other. The int[] is originally in ascending order
        // this generates the next array in lexicographic order descending

        // locate the last occurrence where next[k] < next[k+1]
        int gt = -1;
        for ( int idx = 0; idx < next.length-1; idx++) {
            if ( next[idx] < next[idx+1] ) {
                gt = idx;
            }
        }

        if ( gt == -1 ) {
            return null;
        }

        int largestLessThan = gt+1;
        for ( int idx = 1 + largestLessThan; idx < next.length; idx++) {
            if ( next[gt] < next[idx] ) {
                largestLessThan = idx;
            }
        }

        int val = next[gt];
        next[gt] = next[largestLessThan];
        next[largestLessThan] = val;

        // reverse the tail of the array
        int[] newTail = new int[next.length-gt-1];
        int ctr = 0;
        for ( int idx = next.length-1; idx > gt; idx-- ) {
            newTail[ctr++] = next[idx];
        }

        for ( int idx = 0; idx < newTail.length; idx++) {
            next[gt+idx+1] = newTail[idx];
        }

        return next;
    }

    // man. All this to test dirichlet.

    private double[] unwrap(List<Double> stuff) {
        double[] unwrapped = new double[stuff.size()];
        int idx = 0;
        for ( Double d : stuff ) {
            unwrapped[idx++] = d == null ? 0.0 : d;
        }

        return unwrapped;
    }

    @Test
    public void testNextPermutation() {
        int[] arr = new int[]{1,2,3,4};
        int[][] gens = new int[][] {
                new int[]{1,2,3,4},
                new int[]{1,2,4,3},
                new int[]{1,3,2,4},
                new int[]{1,3,4,2},
                new int[]{1,4,2,3},
                new int[]{1,4,3,2},
                new int[]{2,1,3,4},
                new int[]{2,1,4,3},
                new int[]{2,3,1,4},
                new int[]{2,3,4,1},
                new int[]{2,4,1,3},
                new int[]{2,4,3,1},
                new int[]{3,1,2,4},
                new int[]{3,1,4,2},
                new int[]{3,2,1,4},
                new int[]{3,2,4,1},
                new int[]{3,4,1,2},
                new int[]{3,4,2,1},
                new int[]{4,1,2,3},
                new int[]{4,1,3,2},
                new int[]{4,2,1,3},
                new int[]{4,2,3,1},
                new int[]{4,3,1,2},
                new int[]{4,3,2,1} };
        for ( int gen = 0; gen < gens.length; gen ++ ) {
            for ( int idx = 0; idx < 3; idx++ ) {
                Assert.assertEquals(arr[idx],gens[gen][idx],
                        String.format("Error at generation %d, expected %s, observed %s",gen,Arrays.toString(gens[gen]),Arrays.toString(arr)));
            }
            arr = nextPermutation(arr);
        }
    }

    private double[] addEpsilon(double[] counts) {
        double[] d = new double[counts.length];
        for ( int i = 0; i < counts.length; i ++ ) {
            d[i] = counts[i] + 1e-3;
        }
        return d;
    }

    @Test
    public void testDirichletMultinomial() {
        List<double[]> testAlleles = Arrays.asList(
                new double[]{80,240},
                new double[]{1,10000},
                new double[]{0,500},
                new double[]{5140,20480},
                new double[]{5000,800,200},
                new double[]{6,3,1000},
                new double[]{100,400,300,800},
                new double[]{8000,100,20,80,2},
                new double[]{90,20000,400,20,4,1280,720,1}
        );

        int[] numAlleleSampled = new int[]{2,5,10,20,25};
        for ( double[] alleles : testAlleles ) {
            for ( int count : numAlleleSampled ) {
                // test that everything sums to one. Generate all multinomial draws
                List<Double> likelihoods = new ArrayList<>(100000);
                NextCounts generator = new NextCounts(alleles.length,count);
                double maxLog = Double.MIN_VALUE;
                //List<String> countLog = new ArrayList<String>(200);
                while ( generator.hasNext() ) {
                    int[] thisCount = generator.next();
                    //countLog.add(Arrays.toString(thisCount));
                    Double likelihood = MathUtils.dirichletMultinomial(addEpsilon(alleles),thisCount);
                    Assert.assertTrue(! Double.isNaN(likelihood) && ! Double.isInfinite(likelihood),
                            String.format("Likelihood for counts %s and nAlleles %d was %s",
                                    Arrays.toString(thisCount),alleles.length,Double.toString(likelihood)));
                    if ( likelihood > maxLog )
                        maxLog = likelihood;
                    likelihoods.add(likelihood);
                }
                System.out.printf("%d likelihoods and max is (probability) %e\n",likelihoods.size(),Math.pow(10,maxLog));
                Assert.assertEquals(MathUtils.sumLog10(unwrap(likelihoods)),1.0,1e-7,
                        String.format("Counts %d and alleles %d have nLikelihoods %d. \n Counts: %s",
                                count,alleles.length,likelihoods.size(), "NODEBUG"/*,countLog*/));
            }
        }
    }


    private static final List<double[]> testArrays = Arrays.asList(
            new double[] {1,2,3,4,5},
            new double[] {0.0},
            new double[] {3, 2, 5, 6},
            new double[] {19,-5, 22, 55, -1000, 2,2,2},
            new double[] {-1,-1,-1,-1,-1},
            new double[] {-1,-2,-3,-10,-1}
    );

    @Test
    public void testApplyToArray() {
        for (final double[] array : testArrays) {
            final double[] copy = Arrays.copyOf(array, array.length);
            final DoubleUnaryOperator func = Math::exp;
            Assert.assertEquals(MathUtils.applyToArray(copy, func), Arrays.stream(copy).map(func).toArray());
            Assert.assertEquals(array, copy);   //make sure original array was not affected
        }
    }

    @Test
    public void testApplyToArrayInPlace() {
        for (final double[] array : testArrays) {
            final double[] copy = Arrays.copyOf(array, array.length);
            final DoubleUnaryOperator func = Math::exp;
            final double[] result = MathUtils.applyToArrayInPlace(copy, func);
            Assert.assertTrue(result == copy);
            Assert.assertEquals(copy, Arrays.stream(array).map(func).toArray());   //make sure original array WAS affected
        }
    }

    @Test
    public void testAllMatchDouble() {
        for (final double[] array : testArrays) {
            final double[] copy = Arrays.copyOf(array, array.length);
            final DoublePredicate pred = x -> x > -0.5;
            Assert.assertEquals(MathUtils.allMatch(copy, pred), Arrays.stream(copy).allMatch(pred));
            Assert.assertEquals(array, copy);   //make sure original array was not affected
        }
    }

    @Test
    public void testAllMatchInt() {
        for (final double[] doubleArray : testArrays) {
            final int[] array = Arrays.stream(doubleArray).mapToInt(x -> (int) Math.round(x)).toArray();
            final IntPredicate pred = x -> x > -1;
            Assert.assertEquals(MathUtils.allMatch(array, pred), Arrays.stream(array).allMatch(pred));
        }
    }

    @DataProvider(name = "distanceSquaredData")
    public Object[][] getDistanceSquaredData() {
        return new Object[][] {
                { new double[]{0.0}, new double[]{0.0}, 0.0},
                { new double[]{2.0, 3.0, 4.0}, new double[]{2.0, 3.0, 4.0}, 0.0},
                { new double[]{-2.0, -3.0, -4.0}, new double[]{-2.0, -3.0, -4.0}, 0.0},
                { new double[]{2.0, 3.0, 4.0}, new double[]{3.0, 4.0, 5.0}, 3.0},
        };
    }

    @Test(dataProvider = "distanceSquaredData")
    public void testDistanceSquared(final double[] v1, final double[] v2, final double expected) {
        Assert.assertEquals(MathUtils.distanceSquared(v1, v2), expected);
    }

    @DataProvider(name = "distanceSquaredNegativeData")
    public Object[][] getDistanceSquaredNegativeData() {
        return new Object[][] {
                { new double[]{2.0, 3.0, 4.0}, new double[]{3.0, 4.0, }, Double.NaN},
        };
    }

    @Test(dataProvider = "distanceSquaredNegativeData", expectedExceptions = ArrayIndexOutOfBoundsException.class)
    public void testDistanceSquaredDifferentLengths(final double[] v1, final double[] v2, final double expected) {
        MathUtils.distanceSquared(v1, v2);
        Assert.fail("Vectors with unequal lengths should throw");
    }

    @Test
    public void naturalLogSumNaturalLogTest() {
        final double[] realSpaceValues = {0.1, 0.2, 1.4, 5.9};
        final double[] logSpaceValues= Arrays.stream(realSpaceValues).map(Math::log).toArray();

        final double expected = Math.log(Arrays.stream(realSpaceValues).sum());
        final double actual = NaturalLogUtils.logSumExp(logSpaceValues);
        Assert.assertEquals(actual, expected, 1e-10);
    }

    @Test
    public void testRandomSelectFlatProbability() {
        final RandomGenerator rg = RandomGeneratorFactory.createRandomGenerator(new Random(13));
        final int NUM_SAMPLES = 1000;
        final List<Integer> choices = Arrays.asList(0,1,2);
        final List<Integer> result = IntStream.range(0, NUM_SAMPLES)
                .map(n -> MathUtils.randomSelect(choices, j -> 1.0 / choices.size(), rg))
                .boxed()
                .collect(Collectors.toList());
        Assert.assertEquals(result.stream().filter(n -> n == 0).count(), NUM_SAMPLES / choices.size(), 50);
    }

    @Test
    public void testRandomSelect() {
        final RandomGenerator rg = RandomGeneratorFactory.createRandomGenerator(new Random(13));
        final int NUM_SAMPLES = 1000;
        final List<Integer> choices = Arrays.asList(-1,0,1);
        final List<Integer> result = IntStream.range(0, NUM_SAMPLES)
                .map(n -> MathUtils.randomSelect(choices, j -> j*j/2.0, rg))
                .boxed()
                .collect(Collectors.toList());
        Assert.assertEquals(result.stream().filter(n -> n==0).count(), 0);
        Assert.assertEquals(result.stream().filter(n -> n == 1).count(), NUM_SAMPLES / 2, 50);
    }

    @Test
    public void testSecondSmallestMinusSmallest() {
        Assert.assertEquals(MathUtils.secondSmallestMinusSmallest(new int[0], -1), -1);
        Assert.assertEquals(MathUtils.secondSmallestMinusSmallest(new int[] { 10 }, 3), 3);
        Assert.assertEquals(MathUtils.secondSmallestMinusSmallest(new int[] { -10, -23, 3}, -1), 13);
        Assert.assertEquals(MathUtils.secondSmallestMinusSmallest(new int[] { -10, -23, 3}, -1), 13);
    }

    @Test
    public void testSumArrayFunction(){
        final int min = 0;
        final int max = 5;

        // use 1 + 2 + ... + n = n(n+1)/2
        // [0,1,2] + [1,2,3] + ... + [4,5,6] = [0 + 1 + ... + 4, 1 + ... + 5, 2 + ... + 6] = [10, 15, 21 - 1] = [10, 15, 20]
        final double[] array1 = {10, 15, 20};
        ArrayAsserts.assertArrayEquals(MathUtils.sumArrayFunction(min, max, i -> new double[]{(double) i, i + 1.0, i + 2.0}), array1, 1e-8);
    }

    @Test
    public void testMaxElementIndex(){
        final int[] array1 = new int[]{ 1, 2, 3, 4, 5, 6};
        final int expected1 = 5;
        final int[] array2 = new int[]{ 6, 5, 4, 3, 2, 1};
        final int expected2 = 0;
        final int[] array3 = new int[]{ 3, 6, 2, 1, 4, 5};
        final int expected3 = 1;
        final int[] array4 = new int[]{ 3, 3, 3, 5, 5, 5};
        final int expected4 = 3;

        Assert.assertEquals(MathUtils.maxElementIndex(array1), expected1);
        Assert.assertEquals(MathUtils.maxElementIndex(array2), expected2);
        Assert.assertEquals(MathUtils.maxElementIndex(array3), expected3);
        Assert.assertEquals(MathUtils.maxElementIndex(array4), expected4);
    }

    @DataProvider
    public Object[][] getThrowingCases(){
        return new Object[][]{
                {(long)(Integer.MAX_VALUE) + 1L},
                {(long)(Integer.MIN_VALUE) - 1L},
                {(long)(Integer.MAX_VALUE) * 100L},
        };

    }

    @Test(dataProvider = "getThrowingCases", expectedExceptions = GATKException.class)
    public void testToIntExactOrThrowThrowingCases(long value){
        MathUtils.toIntExactOrThrow(value, () -> new GATKException("can't convert"));
    }

    @DataProvider
    public Object[][] getGoodCases(){
        return new Object[][]{
                {(long)Integer.MAX_VALUE},
                {10000L},
                {10L},
                {0L},
                {-10000L},
                {(long)Integer.MIN_VALUE}
        };
    }

    @Test(dataProvider = "getGoodCases")
    public void testToIntExactOrThrowGoodCases(long value){
        Assert.assertEquals(MathUtils.toIntExactOrThrow(value,() -> new GATKException("")), (int)value);
    }

    @Test
    public void testFastBernoulliEntropy() {
        final int N = 100;
        IntStream.range(0, N + 1).forEach(n -> {
            final double p = (double) n / N;
            final double computed = MathUtils.fastBernoulliEntropy(p);
            if (n == 0 || n == N) {
                Assert.assertEquals(computed, 0, 1.0e-8);
            } else {
                final double exact = -(p * Math.log(p) + (1-p)*Math.log(1-p));
                Assert.assertEquals(computed, exact, 0.01);
            }
        });
    }

    // function, min, max, precision, solution
    @DataProvider
    public Object[][] binarySearchCases(){
        return new Object[][]{
                {(DoubleUnaryOperator) x -> x, -1.0, 1.0, 0.001, OptionalDouble.of(0.0)},
                {(DoubleUnaryOperator) x -> x, 1.0, 2.0, 0.001, OptionalDouble.empty()},
                {(DoubleUnaryOperator) x -> x*x - 2, 0.0, 2.0, 0.001, OptionalDouble.of(Math.sqrt(2.0))},
                {(DoubleUnaryOperator) x -> x*x + 2, 0.0, 2.0, 0.001, OptionalDouble.empty()},
                {(DoubleUnaryOperator) x -> Math.pow((x-4),3), 0, 10, 0.001, OptionalDouble.of(4.0)}
        };
    }

    @Test(dataProvider = "binarySearchCases")
    public void testBinarySearchFindZero(final DoubleUnaryOperator func, final double lower,
                                                      final double upper, final double precision, final OptionalDouble expected) {
        final OptionalDouble result = MathUtils.binarySearchFindZero(func, lower, upper, precision);
        Assert.assertEquals(expected.isPresent(), result.isPresent());
        expected.ifPresent(exp -> Assert.assertEquals(exp, result.getAsDouble(), precision));
    }
}
