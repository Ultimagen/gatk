package org.ultimagen.reads;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.ultimagen.flowBasedRead.read.FlowBasedHaplotype;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class FlowBasedHaplotypeUnitTest extends GATKBaseTest {

    @DataProvider(name="haplotypeGenerator")
    public Object[][] haplotypeTestDataGenerator(){
        String [] haplotypeSeqs = {"ATCGCAGGGAATTGTCCCCATGAAACTAAG",
                                "TGGGCTACCCCGTATATTTCGATTGCATTA",
                                "CCGCCTATTCGCTCTATCGCATCAAATCAA",
                                "GACGGCCTAGCTGCTCGTAGGCATCCTATA",
                                "ACTAACCGCATTTAACGCTCACGCATAAAG",
                                "TGAGTTTTCCACGACGTATTTCAGCTAAGA",
                                "AACTTTCACGTCACACAAGATTCCAGGTAC",
                                "GCACCGTCGTTTCGCCAATAAAATCGACTA",
                                "ATTTCCGCCCCTTAGACATTTGTATAACAT",
                                "GGACTTCGAAATTTTACAGATCATCGCTAC"};
        int[][] expectedFlow =  { {0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 3, 0, 2, 0, 0, 2, 0, 0, 1, 1, 0, 4, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 3, 1, 0, 1, 2, 0, 1},
                {1, 0, 0, 3, 0, 0, 1, 0, 1, 1, 4, 1, 1, 1, 0, 0, 1, 1, 0, 0, 3, 0, 1, 1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 1},
                {0, 0, 2, 1, 0, 0, 2, 0, 1, 1, 0, 0, 2, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 3, 0, 0, 1, 0, 1, 0, 0, 2},
                {0, 0, 0, 1, 0, 1, 1, 2, 0, 0, 2, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 2, 0, 1, 1, 0, 0, 1, 1},
                {0, 1, 1, 0, 1, 2, 2, 1, 0, 0, 1, 0, 0, 1, 0, 0, 3, 2, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 3, 0, 1},
                {1, 0, 0, 1, 0, 1, 0, 1, 4, 0, 2, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 3, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 2, 0, 1, 0, 1},
                {0, 2, 1, 0, 3, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 2, 0, 1, 0, 1, 0, 0, 2, 0, 2, 0, 0, 1, 0, 2, 1, 1, 1},
                {0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 1, 3, 0, 1, 1, 0, 0, 2, 0, 0, 2, 0, 0, 1, 4, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1},
                {0, 1, 0, 0, 3, 0, 2, 1, 0, 0, 4, 0, 2, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 3, 0, 0, 1, 1, 1, 0, 0, 1, 2, 1, 0, 0, 1, 0, 0, 1},
                {0, 0, 0, 2, 0, 1, 1, 0, 2, 0, 1, 1, 0, 3, 0, 0, 4, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1}};


        int [] trimLeft = {0,1,2,3,4,5,6,7,8,9};


        final List<Object[]> tests = new LinkedList<>();
        for (int i = 0; i < haplotypeSeqs.length; i++){
            tests.add( new Object[]{ haplotypeSeqs[i], expectedFlow[i], trimLeft[i], trimLeft[i]});

        }

        return tests.toArray(new Object[][]{});

    }

    @Test(dataProvider = "haplotypeGenerator")
    public void testFlowBasedHaplotype(String inputSeq, int[] expectedKey, int trimLeftIgnore, int trimRightIgnore){
        Haplotype hap = new Haplotype(inputSeq.getBytes());
        FlowBasedHaplotype fbh = new FlowBasedHaplotype(hap, "TACG");
        Assert.assertEquals(fbh.getKey(), expectedKey);


    }
    @Test(dataProvider = "haplotypeGenerator")
    public void testFindLeftClipping(String inputSeq, int[] expectedKey, int trimLeft, int trimRightIgnore) {
        Haplotype hap = new Haplotype(inputSeq.getBytes());
        Haplotype hapTrimmed = new Haplotype(inputSeq.substring(trimLeft).getBytes());
        FlowBasedHaplotype fbh = new FlowBasedHaplotype(hap, "TACG");
        int[] leftClip = fbh.findLeftClipping(trimLeft);
        int[] trimmedKey= Arrays.copyOfRange(fbh.getKey(), leftClip[0], fbh.getKeyLength());
        trimmedKey[0]-=leftClip[1];
        FlowBasedHaplotype fbhTrimmed = new FlowBasedHaplotype(hapTrimmed, "TACG");
        if (trimLeft>0) {
            expectedKey = Arrays.copyOfRange(fbhTrimmed.getKey(), findFirstNonzero(fbhTrimmed.getKey()), fbhTrimmed.getKeyLength());
        } else {
            expectedKey = fbhTrimmed.getKey();
        }
        Assert.assertEquals(trimmedKey, expectedKey);

    }


    @Test(dataProvider = "haplotypeGenerator")
    public void testFindRightClipping(String inputSeq, int[] expectedKey, int trimLeftIgnore, int trimRight) {
        Haplotype hap = new Haplotype(inputSeq.getBytes());
        Haplotype hapTrimmed = new Haplotype(inputSeq.substring(0, inputSeq.length() - trimRight).getBytes());
        FlowBasedHaplotype fbh = new FlowBasedHaplotype(hap, "TACG");
        int[] rightClip = fbh.findRightClipping(trimRight);
        int[] trimmedKey= Arrays.copyOfRange(fbh.getKey(), 0, fbh.getKeyLength()-rightClip[0]);
        trimmedKey[trimmedKey.length-1]-=rightClip[1];
        FlowBasedHaplotype fbhTrimmed = new FlowBasedHaplotype(hapTrimmed, "TACG");
        if (trimRight>0) {
            expectedKey = Arrays.copyOfRange(fbhTrimmed.getKey(), 0, findLastNonzero(fbhTrimmed.getKey())+1);
        } else {
            expectedKey = fbhTrimmed.getKey();
        }
        Assert.assertEquals(trimmedKey, expectedKey);

    }

    @DataProvider(name="haplotypeModificationProvider")
    public Object[][] haplotypeModifier(){
        boolean[] answers = {true, true, true, false, false};
        String sourceHaplotype = "ATCGCAGGGAATTGTCCCCATGAAACTAAG";
        String[] modifiedHaplotypes = {"ATCGCAGGGAATTGTCCCCATGAAACTAAG",
                "ATCGCAGGGGAATTGTCCCCATGAAACTAAG",
                "ATCGCAGGGATTGTCCCCATGAAACTAAG",
                "ATCGCAGGGTTGTCCCCATGAAACTAAG",
                "ATCGCAGGGAATAGTCCCCATGAAACTAAG"};


        final List<Object[]> tests = new LinkedList<>();
        for (int i = 0; i < modifiedHaplotypes.length; i++){
            tests.add( new Object[]{ new Haplotype(sourceHaplotype.getBytes()),
                    new Haplotype(modifiedHaplotypes[i].getBytes()), answers[i]});

        }

        return tests.toArray(new Object[][]{});

    }


    @Test(dataProvider = "haplotypeModificationProvider")
    public void testEqualUpToHmerChange(Haplotype hap1, Haplotype hap2, boolean answer) {
        FlowBasedHaplotype fbh1 = new FlowBasedHaplotype(hap1, "TACG");
        FlowBasedHaplotype fbh2 = new FlowBasedHaplotype(hap2, "TACG");
        Assert.assertEquals(fbh1.equalUpToHmerChange(fbh2),answer);
    }

    private int findFirstNonzero(int [] key){
        for (int i = 0; i < key.length; i++){
            if (key[i]>0){
                return i;
            }
        }
        return key.length;
    }

    private int findLastNonzero(int [] key){
        for (int i = key.length-1; i >= 0; i--){
            if (key[i]>0){
                return i;
            }
        }
        return -1;
    }

}