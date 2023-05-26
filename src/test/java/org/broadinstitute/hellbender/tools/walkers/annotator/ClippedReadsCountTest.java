package org.broadinstitute.hellbender.tools.walkers.annotator;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.testng.Assert.*;

public class ClippedReadsCountTest extends GATKBaseTest {
    private static final Allele REF = Allele.create("A", true);
    private static final Allele ALT = Allele.create("C");
    private static final List<Allele> ALLELES = Arrays.asList(REF, ALT);
    private static final String SAMPLE = "sample1";

    private GATKRead makeRead() {
        return ArtificialAnnotationUtils.makeRead(30, 50);
    }
    private GATKRead makeSoftClippedRead(final boolean leftClip){ return ArtificialAnnotationUtils.makeSoftClippedRead(10, 50, leftClip);}

    @Test
    public void testDescription(){
        String[] constants = {GATKVCFConstants.SOFT_CLIP_LEFT_COUNT_KEY, GATKVCFConstants.SOFT_CLIP_RIGHT_COUNT_KEY};
        VCFFormatHeaderLine[] hlines    = {GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.SOFT_CLIP_LEFT_COUNT_KEY),
                              GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.SOFT_CLIP_RIGHT_COUNT_KEY)};
        Assert.assertEquals(new ClippedReadsCount().getKeyNames(), new ArrayList<>(Arrays.asList(constants)));
        Assert.assertEquals(new ClippedReadsCount().getDescriptions(), new ArrayList<>(Arrays.asList(hlines)));
    }

    @DataProvider(name = "SoftClippingData")
    public Object[][][] readDepthData() {
        return new Object[][][]{{
            {20,0}, {0,17}},
            {{0,0}, {0,0}}
        };
    }

    @Test(dataProvider = "SoftClippingData")
    public void testUsingReads(final Object[] refDepth, final Object[] altDepth){
        final int[] expectedLeftCount = {(int)refDepth[0], (int)altDepth[0]};
        final int[] expectedRightCount = {(int)refDepth[1], (int)altDepth[1]};


        final int dpDepth = 50; //Note: using a different value on purpose so that we can check that reads are preferred over DP
        final Genotype gAC = new GenotypeBuilder(SAMPLE, ALLELES).DP(dpDepth).make();

        final double log10PError = -5;

        final List<GATKRead> refReads = IntStream.range(0, dpDepth - expectedLeftCount[0] - expectedRightCount[0]).mapToObj(i -> makeRead()).collect(Collectors.toList());
        final List<GATKRead> refClippedReads = IntStream.range(0, expectedLeftCount[0]).mapToObj(i -> makeSoftClippedRead(true)).collect(Collectors.toList());
        final List<GATKRead> altClippedReads = IntStream.range(0, expectedRightCount[1]).mapToObj(i -> makeSoftClippedRead(false)).collect(Collectors.toList());
        refReads.addAll(refClippedReads);
        final AlleleLikelihoods<GATKRead, Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(SAMPLE, refReads, altClippedReads, -100.0, -100.0, REF, ALT);

        final VariantContext vc = new VariantContextBuilder("test", "20", 10, 10, ALLELES).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();

        final GenotypeBuilder gb = new GenotypeBuilder(gAC);
        new ClippedReadsCount().annotate(new ReferenceContext(new ReferenceFileSource(Path.of(b37_reference_20_21)),
                new SimpleInterval("20", 10, 11)), vc, gAC, gb, likelihoods);
        final int[] ad = gb.make().getAD();
        Assert.assertEquals(ad, new int[]{refReads.size(), altClippedReads.size()});

        //now test a no-op
        final GenotypeBuilder gb1 = new GenotypeBuilder(gAC);
        new ClippedReadsCount().annotate(null, vc, null, gb1, likelihoods);  //null genotype
        Assert.assertFalse(gb1.make().hasAD());
    }
}