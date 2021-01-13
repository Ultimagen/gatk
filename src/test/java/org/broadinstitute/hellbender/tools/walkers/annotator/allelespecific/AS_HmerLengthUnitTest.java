package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Map;

public class AS_HmerLengthUnitTest extends GATKBaseTest {

    private final ReferenceDataSource refSource = ReferenceDataSource.of(new File("src/test/resources/hg19micro.fasta").toPath());

    // REF:  T A A - C
    // ALT1: T A A A C HmerLength: 3
    // ALT2: T A - - C HmerLength: 2
    // ALT3: C A A - C HmerLength: 0
    // ALT4: T - - - C HmerLength: 2
    @Test
    public void testHmerAlleles() {
        VariantContext vc = new VariantContextBuilder().chr("1").start(10025).stop(10027).alleles("TAA", "TAAA", "TA", "CAA", "T").make();
        ReferenceContext ref = new ReferenceContext(refSource, new SimpleInterval("1", 10025, 10125));
        AS_HmerLength annotation = new AS_HmerLength();
        Map<String, Object> output = annotation.annotate(ref, vc, null);
        Assert.assertEquals(output.get(GATKVCFConstants.AS_HMER_LENGTH_KEY), "3,2,0,2");
    }

    // REF:  T A A C C C T A A C C C
    // ALT1: T - - - - - - A A C C C HmerLength: 0
    // ALT2: C T A C C C T A A C C C HmerLenth: 0 //MNP should be length 0?
    @Test
    public void testNonHmerAlleles() {
        VariantContext vc = new VariantContextBuilder().chr("1").start(10025).stop(10031).alleles("TAACCCT", "T", "CTACCCT").make();
        ReferenceContext ref = new ReferenceContext(refSource, new SimpleInterval("1", 10025, 10031));
        AS_HmerLength annotation = new AS_HmerLength();
        Map<String, Object> output = annotation.annotate(ref, vc, null);
        Assert.assertEquals(output.get(GATKVCFConstants.AS_HMER_LENGTH_KEY), "0,0");
    }

    // REF:  C C C G C C C
    // ALT1: C C C C C C C HmerLength: 0 since this is a SNP
    @Test
    public void testUnhandledHmers() {
        VariantContext vc = new VariantContextBuilder().chr("2").start(10803).stop(10803).alleles("G", "C").make();
        ReferenceContext ref = new ReferenceContext(refSource, new SimpleInterval("2", 10803, 10803));
        AS_HmerLength annotation = new AS_HmerLength();
        Map<String, Object> output = annotation.annotate(ref, vc, null);
        Assert.assertEquals(output.get(GATKVCFConstants.AS_HMER_LENGTH_KEY), "0");
    }

    // REF:  G G G G T C C C
    // ALT1: G G G G C C C C
    // ALT2: G G G G G C C C
    @Test
    public void testAllSnps() {
        VariantContext vc = new VariantContextBuilder().chr("2").start(10679).stop(10679).alleles("T", "C", "G").make();
        ReferenceContext ref = new ReferenceContext(refSource, new SimpleInterval("2", 10679, 10679));
        AS_HmerLength annotation = new AS_HmerLength();
        Map<String, Object> output = annotation.annotate(ref, vc, null);
        Assert.assertEquals(output.get(GATKVCFConstants.AS_HMER_LENGTH_KEY), "0,0");
    }

    // REF:  G G - G G T C C C
    // ALT1: G G - - G T C C C //not left aligned deletion
    // ALT2: G G G G G T C C C //not left aligned since there are G's to the left of this G insertion
    // ALT3: G G T G G T C C C // left aligned because it is anchored in non h-mer bases
    // ALT4: G G - - - - - C C //left aligned because it is anchored in non h-mer bases
    @Test
    public void testNonLeftAlignedIndel() {
        VariantContext notLeftAlignedVc = new VariantContextBuilder().chr("2").start(10676).stop(10677).alleles("GG", "G", "GGG").make();
        ReferenceContext ref = new ReferenceContext(refSource, new SimpleInterval("2", 10676, 10680));
        AS_HmerLength annotation = new AS_HmerLength();
        Map<String, Object> notLeftAlignedOutput = annotation.annotate(ref, notLeftAlignedVc, null);
        Assert.assertEquals(notLeftAlignedOutput.get(GATKVCFConstants.AS_HMER_LENGTH_KEY), "4,5");
        VariantContext leftAlignedVc = new VariantContextBuilder().chr("2").start(10676).stop(10680).alleles("GGGTC", "GTGGTC", "G").make();
        Map<String, Object> output = annotation.annotate(ref, leftAlignedVc, null);
        Assert.assertEquals(output.get(GATKVCFConstants.AS_HMER_LENGTH_KEY), "1,0");
    }

    // REF:  G G G G T - C C C
    // ALT1: G G G G T T C C C //not left aligned because it should be:

    // REF:  G G G G - T C C C
    // ALT1: G G G G T T C C C
    @Test
    public void testNonLeftAlignedInsertion() {
        VariantContext nonLeftAlignedVc = new VariantContextBuilder().chr("2").start(10679).stop(10679).alleles("T", "TT").make();
        ReferenceContext ref = new ReferenceContext(refSource, new SimpleInterval("2", 10679, 10680));
        AS_HmerLength annotation = new AS_HmerLength();
        Map<String, Object> output = annotation.annotate(ref, nonLeftAlignedVc, null);
        Assert.assertEquals(output.get(GATKVCFConstants.AS_HMER_LENGTH_KEY), "2");
    }

    // REF:  G G G G - T C C C
    // ALT1: G G G G A T C C C //hmer length: 1
    @Test
    public void testRightAnchoredInsertion() {
        VariantContext vc = new VariantContextBuilder().chr("2").start(10678).stop(10678).alleles("G", "GA").make();
        ReferenceContext ref = new ReferenceContext(refSource, new SimpleInterval("2", 10678, 10679));
        AS_HmerLength annotation = new AS_HmerLength();
        Map<String, Object> output = annotation.annotate(ref, vc, null);
        Assert.assertEquals(output.get(GATKVCFConstants.AS_HMER_LENGTH_KEY), "1");
    }

    // REF:        T C C C A
    // REF allele: T C C C A   //untrimmed
    // ALT allele: T C - - A //hmer length: 3
    @Test
    public void testUntrimmed() {
        VariantContext vc = new VariantContextBuilder().chr("2").start(10679).stop(10682).alleles("TCCC", "TC").make();
        ReferenceContext ref = new ReferenceContext(refSource, new SimpleInterval("2", 10679, 10682));
        AS_HmerLength annotation = new AS_HmerLength();
        Map<String, Object> output = annotation. annotate(ref, vc, null);
        Assert.assertEquals(output.get(GATKVCFConstants.AS_HMER_LENGTH_KEY), "3");
    }
}
