package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.mockito.Mockito;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class FlowBasedPileupElementTest extends GATKBaseTest {

    @DataProvider(name = "flowBasedPileupData")
    public Object[][] getFlowBasedPileupData() {
        return new Object[][] {
                { 0, 0, (byte)'A', 3, new double[] {0.9, 0.1} },
                { 1, 2, (byte)'C', 5, new double[] {0.8, 0.2} },
                { 2, 5, (byte)'G', 7, new double[] {0.7, 0.3} }
        };
    }

    /**
     * Helper method to create a real PileupElement for testing
     */
    private PileupElement createPileupElement(GATKRead read, int offset) {
        // Create a simple CigarElement for testing (1M - one match)
        CigarElement cigarElement = new CigarElement(1, CigarOperator.M);
        return new PileupElement(read, offset, cigarElement, 0, 0);
    }

    @Test(dataProvider = "flowBasedPileupData")
    public void getFlowNumReturnsCorrectValue(int offset, int expectedFlowNum, byte expectedNuc, int expectedHpolCall, double[] expectedProbs) {
        FlowBasedRead mockRead = Mockito.mock(FlowBasedRead.class);
        GATKRead mockGATKRead = Mockito.mock(GATKRead.class);
        PileupElement pileup = createPileupElement(mockGATKRead, offset);

        Mockito.when(mockRead.getBase2Flow(offset)).thenReturn(expectedFlowNum);

        FlowBasedPileupElement element = new FlowBasedPileupElement(mockRead, pileup, false);
        Assert.assertEquals(element.getFlowNum(), expectedFlowNum);
    }

    @Test(dataProvider = "flowBasedPileupData")
    public void getFlowNucReturnsCorrectNucleotide(int offset, int flowNum, byte expectedNuc, int expectedHpolCall, double[] expectedProbs) {
        FlowBasedRead mockRead = Mockito.mock(FlowBasedRead.class);
        GATKRead mockGATKRead = Mockito.mock(GATKRead.class);
        PileupElement pileup = createPileupElement(mockGATKRead, offset);

        Mockito.when(mockRead.getBase2Flow(offset)).thenReturn(flowNum);
        Mockito.when(mockRead.getNucForFlow(flowNum)).thenReturn(expectedNuc);

        FlowBasedPileupElement element = new FlowBasedPileupElement(mockRead, pileup, false);
        Assert.assertEquals(element.getFlowNuc(), expectedNuc);
    }

    @Test(dataProvider = "flowBasedPileupData")
    public void getHpolCallReturnsCorrectValue(int offset, int flowNum, byte expectedNuc, int expectedHpolCall, double[] expectedProbs) {
        FlowBasedRead mockRead = Mockito.mock(FlowBasedRead.class);
        GATKRead mockGATKRead = Mockito.mock(GATKRead.class);
        PileupElement pileup = createPileupElement(mockGATKRead, offset);
        int[] keys = {3, 5, 7, 9, 2, 1};

        Mockito.when(mockRead.getBase2Flow(offset)).thenReturn(flowNum);
        Mockito.when(mockRead.getKey()).thenReturn(keys);

        FlowBasedPileupElement element = new FlowBasedPileupElement(mockRead, pileup, false);
        Assert.assertEquals(element.getHpolCall(), keys[flowNum]);
    }

    @Test(dataProvider = "flowBasedPileupData")
    public void getCallProbsReturnsCorrectValues(int offset, int flowNum, byte expectedNuc, int expectedHpolCall, double[] expectedProbs) {
        FlowBasedRead mockRead = Mockito.mock(FlowBasedRead.class);
        GATKRead mockGATKRead = Mockito.mock(GATKRead.class);
        PileupElement pileup = createPileupElement(mockGATKRead, offset);

        Mockito.when(mockRead.getBase2Flow(offset)).thenReturn(flowNum);
        Mockito.when(mockRead.getHmerProbs(flowNum)).thenReturn(expectedProbs);

        FlowBasedPileupElement element = new FlowBasedPileupElement(mockRead, pileup, false);
        Assert.assertEquals(element.getCallProbs(), expectedProbs);
    }


    @Test(expectedExceptions = ArrayIndexOutOfBoundsException.class)
    public void getFlowNumThrowsWithInvalidOffset() {
        FlowBasedRead mockRead = Mockito.mock(FlowBasedRead.class);
        GATKRead mockGATKRead = Mockito.mock(GATKRead.class);
        PileupElement pileup = createPileupElement(mockGATKRead, 10);

        Mockito.when(mockRead.getBase2Flow(10)).thenThrow(new ArrayIndexOutOfBoundsException("Invalid offset"));

        FlowBasedPileupElement element = new FlowBasedPileupElement(mockRead, pileup, false);
        element.getFlowNum();
    }

    @Test
    public void getHpolCallHandlesEmptyKeyArray() {
        FlowBasedRead mockRead = Mockito.mock(FlowBasedRead.class);
        GATKRead mockGATKRead = Mockito.mock(GATKRead.class);
        PileupElement pileup = createPileupElement(mockGATKRead, 0);
        int[] emptyKeys = {};

        Mockito.when(mockRead.getBase2Flow(0)).thenReturn(0);
        Mockito.when(mockRead.getKey()).thenReturn(emptyKeys);

        FlowBasedPileupElement element = new FlowBasedPileupElement(mockRead, pileup, false);

        try {
            element.getHpolCall();
            Assert.fail("Expected ArrayIndexOutOfBoundsException");
        } catch (ArrayIndexOutOfBoundsException e) {
            // Expected exception
        }
    }

    @Test
    public void getCallProbsHandlesEmptyProbArray() {
        FlowBasedRead mockRead = Mockito.mock(FlowBasedRead.class);
        GATKRead mockGATKRead = Mockito.mock(GATKRead.class);
        PileupElement pileup = createPileupElement(mockGATKRead, 0);
        double[] emptyProbs = {};

        Mockito.when(mockRead.getBase2Flow(0)).thenReturn(0);
        Mockito.when(mockRead.getHmerProbs(0)).thenReturn(emptyProbs);

        FlowBasedPileupElement element = new FlowBasedPileupElement(mockRead, pileup, false);
        Assert.assertEquals(element.getCallProbs().length, 0);
    }
}
