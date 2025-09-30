package org.broadinstitute.hellbender.utils.pileup;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
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

    @Test(dataProvider = "flowBasedPileupData")
    public void getFlowNumReturnsCorrectValue(int offset, int expectedFlowNum, byte expectedNuc, int expectedHpolCall, double[] expectedProbs) {
        FlowBasedRead mockRead = Mockito.mock(FlowBasedRead.class);
        PileupElement mockPileup = Mockito.mock(PileupElement.class);

        Mockito.when(mockPileup.getOffset()).thenReturn(offset);
        Mockito.when(mockRead.getBase2Flow(offset)).thenReturn(expectedFlowNum);

        FlowBasedPileupElement element = new FlowBasedPileupElement(mockRead, mockPileup);
        Assert.assertEquals(element.getFlowNum(), expectedFlowNum);
    }

    @Test(dataProvider = "flowBasedPileupData")
    public void getFlowNucReturnsCorrectNucleotide(int offset, int flowNum, byte expectedNuc, int expectedHpolCall, double[] expectedProbs) {
        FlowBasedRead mockRead = Mockito.mock(FlowBasedRead.class);
        PileupElement mockPileup = Mockito.mock(PileupElement.class);

        Mockito.when(mockPileup.getOffset()).thenReturn(offset);
        Mockito.when(mockRead.getBase2Flow(offset)).thenReturn(flowNum);
        Mockito.when(mockRead.getNucForFlow(flowNum)).thenReturn(expectedNuc);

        FlowBasedPileupElement element = new FlowBasedPileupElement(mockRead, mockPileup);
        Assert.assertEquals(element.getFlowNuc(), expectedNuc);
    }

    @Test(dataProvider = "flowBasedPileupData")
    public void getHpolCallReturnsCorrectValue(int offset, int flowNum, byte expectedNuc, int expectedHpolCall, double[] expectedProbs) {
        FlowBasedRead mockRead = Mockito.mock(FlowBasedRead.class);
        PileupElement mockPileup = Mockito.mock(PileupElement.class);
        int[] keys = {3, 5, 7, 9, 2, 1};

        Mockito.when(mockPileup.getOffset()).thenReturn(offset);
        Mockito.when(mockRead.getBase2Flow(offset)).thenReturn(flowNum);
        Mockito.when(mockRead.getKey()).thenReturn(keys);

        FlowBasedPileupElement element = new FlowBasedPileupElement(mockRead, mockPileup);
        Assert.assertEquals(element.getHpolCall(), keys[flowNum]);
    }

    @Test(dataProvider = "flowBasedPileupData")
    public void getCallProbsReturnsCorrectValues(int offset, int flowNum, byte expectedNuc, int expectedHpolCall, double[] expectedProbs) {
        FlowBasedRead mockRead = Mockito.mock(FlowBasedRead.class);
        PileupElement mockPileup = Mockito.mock(PileupElement.class);

        Mockito.when(mockPileup.getOffset()).thenReturn(offset);
        Mockito.when(mockRead.getBase2Flow(offset)).thenReturn(flowNum);
        Mockito.when(mockRead.getHmerProbs(flowNum)).thenReturn(expectedProbs);

        FlowBasedPileupElement element = new FlowBasedPileupElement(mockRead, mockPileup);
        Assert.assertEquals(element.getCallProbs(), expectedProbs);
    }


    @Test(expectedExceptions = ArrayIndexOutOfBoundsException.class)
    public void getFlowNumThrowsWithInvalidOffset() {
        FlowBasedRead mockRead = Mockito.mock(FlowBasedRead.class);
        PileupElement mockPileup = Mockito.mock(PileupElement.class);

        Mockito.when(mockPileup.getOffset()).thenReturn(10);
        Mockito.when(mockRead.getBase2Flow(10)).thenThrow(new ArrayIndexOutOfBoundsException("Invalid offset"));

        FlowBasedPileupElement element = new FlowBasedPileupElement(mockRead, mockPileup);
        element.getFlowNum();
    }

    @Test
    public void getHpolCallHandlesEmptyKeyArray() {
        FlowBasedRead mockRead = Mockito.mock(FlowBasedRead.class);
        PileupElement mockPileup = Mockito.mock(PileupElement.class);
        int[] emptyKeys = {};

        Mockito.when(mockPileup.getOffset()).thenReturn(0);
        Mockito.when(mockRead.getBase2Flow(0)).thenReturn(0);
        Mockito.when(mockRead.getKey()).thenReturn(emptyKeys);

        FlowBasedPileupElement element = new FlowBasedPileupElement(mockRead, mockPileup);

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
        PileupElement mockPileup = Mockito.mock(PileupElement.class);
        double[] emptyProbs = {};

        Mockito.when(mockPileup.getOffset()).thenReturn(0);
        Mockito.when(mockRead.getBase2Flow(0)).thenReturn(0);
        Mockito.when(mockRead.getHmerProbs(0)).thenReturn(emptyProbs);

        FlowBasedPileupElement element = new FlowBasedPileupElement(mockRead, mockPileup);
        Assert.assertEquals(element.getCallProbs().length, 0);
    }
}