package org.broadinstitute.hellbender.tools.walkers.genotyper;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.SimpleAllele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.utils.dragstr.DragstrReferenceAnalyzer;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

public final class FlowBasedGenotypesModelUnitTest {
    @DataProvider(name = "isEligibleHomopolymerIndelData")
    public Object[][] isEligibleHomopolymerIndelData() {
        //repeat AG up to length 100, then 6As, then AG to length 200
        String refHap = "AG".repeat(49) + "GG" + "A".repeat(6) + "GA".repeat(100);
        DragstrReferenceAnalyzer dragstr = DragstrReferenceAnalyzer.of(refHap.getBytes(), 0, refHap.length(),1 );
        return new Object[][]{
                {null, 100, dragstr, 5, false},
                {new VariantContextBuilder("test", "1", 100, 100, Arrays.asList(Allele.create("G", true), Allele.create("A", false))).make(), 100, dragstr, 5, false},
                {new VariantContextBuilder("test", "1", 100, 100, Arrays.asList(Allele.create("G", true), Allele.create("GA", false))).make(), 100, dragstr, 5, true},
                {new VariantContextBuilder("test", "1", 100, 100, Arrays.asList(Allele.create("G", true), Allele.create("GA", false))).make(), 100, dragstr, 10, false},
                {new VariantContextBuilder("test", "1", 99, 99, Arrays.asList(Allele.create("G", true), Allele.create("GG", false))).make(), 100, dragstr, 5, false}
        };
    }

    @Test(dataProvider = "isEligibleHomopolymerIndelData")
    public void isEligibleHomopolymerIndel(VariantContext vc, int loc, DragstrReferenceAnalyzer dragstrs, int hpolIndelThreshold, boolean expected) {
        Assert.assertEquals(FlowBasedGenotypesModel.isEligibleHomopolymerIndel(vc, loc, dragstrs, hpolIndelThreshold), expected);
    }

    @DataProvider(name = "isHmerIndelData")
    public Object[][] isHmerIndelData() {
        return new Object[][]{
                {Allele.create("A", false), (byte) 'A', true},
                {Allele.create("AA", false), (byte) 'A', true},
                {Allele.create("AC", false), (byte) 'A', false},
                {Allele.create("AAA", false), (byte) 'A', true},
                {Allele.create("AAT", false), (byte) 'A', false}
        };
    }

    @Test(dataProvider = "isHmerIndelData")
    public void isHmerIndel(Allele al, byte hmer_base, boolean expected) {
        Assert.assertEquals(FlowBasedGenotypesModel.isHmerIndel(al, hmer_base), expected);
    }

}