package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoodsUnitTester;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Tests {@link GenotypesCache}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class GenotypeLikelihoodCalculatorUnitTest {
    
    @Test(dataProvider = "ploidyAndMaximumAlleleAndReadCountsData")
    public void testLikelihoodCalculation(final int ploidy, final int alleleCount, final int[] readCount) {
        final AlleleLikelihoods<GATKRead, Allele> readLikelihoods = ReadLikelihoodsUnitTester.readLikelihoods(alleleCount, readCount);
        final int sampleCount = readCount.length;
        for (int s = 0; s < sampleCount ; s++) {
            final LikelihoodMatrix<GATKRead, Allele> sampleLikelihoods = readLikelihoods.sampleMatrix(s);
            final GenotypeLikelihoods genotypeLikelihoods = GenotypeLikelihoodCalculator.log10GenotypeLikelihoods(ploidy, sampleLikelihoods);
            final double[] genotypeLikelihoodsDoubles = genotypeLikelihoods.getAsVector();
            for (final GenotypeAlleleCounts gac : GenotypeAlleleCounts.iterable(ploidy, alleleCount)) {
                Assert.assertNotNull(genotypeLikelihoods);
                final double[] readGenotypeLikelihoods = new double[sampleLikelihoods.evidenceCount()];
                for (int r = 0; r < sampleLikelihoods.evidenceCount(); r++) {
                    final double[] compoments = new double[gac.distinctAlleleCount()];
                    for (int ar = 0; ar < gac.distinctAlleleCount(); ar++) {
                        final int a = gac.alleleIndexAt(ar);
                        final int aCount = gac.alleleCountAt(ar);
                        final double readLk = sampleLikelihoods.get(a, r);
                        compoments[ar] = readLk + Math.log10(aCount);
                    }
                    readGenotypeLikelihoods[r] = MathUtils.approximateLog10SumLog10(compoments) - Math.log10(ploidy);
                }
                final double genotypeLikelihood = MathUtils.sum(readGenotypeLikelihoods);
                Assert.assertEquals(genotypeLikelihoodsDoubles[gac.index()], genotypeLikelihood, 0.0001 * Math.abs(genotypeLikelihood));
            }
        }
    }

    @Test(dataProvider = "ploidyAndMaximumAlleleAndReadCountsData")
    public void testLikelihoodCalculationWithReferenceBias(final int ploidy, final int alleleCount, final int[] readCount) {
        final AlleleLikelihoods<GATKRead, Allele> readLikelihoods = ReadLikelihoodsUnitTester.readLikelihoods(alleleCount, readCount);
        final int sampleCount = readCount.length;
        final double referenceBias = 0.7;  // Test with a specific reference bias value

        for (int s = 0; s < sampleCount; s++) {
            final LikelihoodMatrix<GATKRead, Allele> sampleLikelihoods = readLikelihoods.sampleMatrix(s);
            final GenotypeLikelihoods genotypeLikelihoods = GenotypeLikelihoodCalculator.log10GenotypeLikelihoods(ploidy, sampleLikelihoods, referenceBias);
            final double[] genotypeLikelihoodsDoubles = genotypeLikelihoods.getAsVector();

            // Calculate expected likelihoods with reference bias
            final double scaledRefBias = referenceBias * ploidy;
            final double altBias = ploidy > 1 ? (ploidy - scaledRefBias) / (ploidy - 1) : 1 - scaledRefBias;
            final double log10RefBias = Math.log10(scaledRefBias);
            final double log10AltBias = Math.log10(altBias);

            for (final GenotypeAlleleCounts gac : GenotypeAlleleCounts.iterable(ploidy, alleleCount)) {
                Assert.assertNotNull(genotypeLikelihoods);
                final double[] readGenotypeLikelihoods = new double[sampleLikelihoods.evidenceCount()];

                boolean containsRefAllele = gac.alleleCountFor(0) > 0;
                boolean containsAltAllele = gac.alleleCountFor(0) < ploidy;

                for (int r = 0; r < sampleLikelihoods.evidenceCount(); r++) {
                    final double[] components = new double[gac.distinctAlleleCount()];
                    for (int ar = 0; ar < gac.distinctAlleleCount(); ar++) {
                        final int a = gac.alleleIndexAt(ar);
                        final int aCount = gac.alleleCountAt(ar);
                        final double readLk = sampleLikelihoods.get(a, r);

                        // Apply bias only to heterozygous genotypes (containing both ref and alt)
                        if (containsRefAllele && containsAltAllele) {
                            final double bias = (a == 0) ? log10RefBias : log10AltBias;
                            components[ar] = readLk + Math.log10(aCount) + bias;
                        } else {
                            components[ar] = readLk + Math.log10(aCount);
                        }
                    }
                    readGenotypeLikelihoods[r] = MathUtils.approximateLog10SumLog10(components) - Math.log10(ploidy);
                }
                final double genotypeLikelihood = MathUtils.sum(readGenotypeLikelihoods);
                Assert.assertEquals(genotypeLikelihoodsDoubles[gac.index()], genotypeLikelihood, 0.0001 * Math.abs(genotypeLikelihood),
                        "Mismatch for genotype " + gac + " with reference bias " + referenceBias);
            }
        }
    }

    @Test
    public void testReferenceBiasAffectsHeterozygousGenotypes() {
        // Test that reference bias changes het genotype likelihoods but not hom genotypes
        final int ploidy = 2;
        final int alleleCount = 2;
        final int[] readCount = {10};

        final AlleleLikelihoods<GATKRead, Allele> readLikelihoods = ReadLikelihoodsUnitTester.readLikelihoods(alleleCount, readCount);
        final LikelihoodMatrix<GATKRead, Allele> sampleLikelihoods = readLikelihoods.sampleMatrix(0);

        // Calculate with default bias (1/ploidy = 0.5)
        final GenotypeLikelihoods defaultLikelihoods = GenotypeLikelihoodCalculator.log10GenotypeLikelihoods(ploidy, sampleLikelihoods);
        final double[] defaultArray = defaultLikelihoods.getAsVector();

        // Calculate with higher reference bias
        final double highRefBias = 0.8;
        final GenotypeLikelihoods biasedLikelihoods = GenotypeLikelihoodCalculator.log10GenotypeLikelihoods(ploidy, sampleLikelihoods, highRefBias);
        final double[] biasedArray = biasedLikelihoods.getAsVector();

        // For diploid with 2 alleles: genotype order is AA(0), AB(1), BB(2)
        // Homozygous genotypes (AA, BB) should be the same
        Assert.assertEquals(defaultArray[0], biasedArray[0], 1e-10, "Hom ref (AA) should not change with reference bias");
        Assert.assertEquals(defaultArray[2], biasedArray[2], 1e-10, "Hom alt (BB) should not change with reference bias");

        // Heterozygous genotype (AB) should be different
        Assert.assertNotEquals(defaultArray[1], biasedArray[1], "Het (AB) should change with reference bias");
    }


    private static final int[] MAXIMUM_ALLELE = { 1, 2, 5, 6};

    private static final int[] PLOIDY = { 1, 2, 3, 20 };

    private static final int[][] READ_COUNTS = {
            { 10 , 100, 50 },
            { 0, 100, 10, 1 , 50 },
            { 1, 2, 3, 4, 20 },
            { 10, 0 },
    };

    @DataProvider(name="ploidyAndMaximumAlleleAndReadCountsData")
    public Object[][] ploidyAndMaximumAlleleAndReadCountsData() {
        final Object[][] result = new Object[PLOIDY.length * MAXIMUM_ALLELE.length * READ_COUNTS.length][];
        int index = 0;
        for (final int i : PLOIDY)
            for (final int j : MAXIMUM_ALLELE)
                for (final int[] k : READ_COUNTS)
                result[index++] = new Object[] { i, j, k };
        return result;
    }
}
