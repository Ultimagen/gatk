package org.ultimagen.flowBasedRead.read;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.ultimagen.variantRecalling.TrimmedReadsReaderUnitTest;

import java.util.ArrayList;
import java.util.List;

public class FlowBasedKeyCodecUnitTest extends GATKBaseTest {

    @DataProvider(name = "testData")
    public Object[][] getTestData() {

        final Object[][]        testData = {

                // trivial cases
                { "T", "TGCA", "1", null },
                { "TT", "TGCA", "2", null },
                { "TGCA", "TGCA", "1,1,1,1", null },
                { "TA", "TGCA", "1,0,0,1", null },
                { "TTAATG", "TGCA", "2,0,0,2,1,1", null },

                // clipping
                { "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                          , "TGCA", "130", "127" },

                // N processing
                { "TNTA", "TGCA", "3,0,0,1", null},
                { "TTNA", "TGCA", "3,0,0,1", null},
                { "TTAN", "TGCA", "2,0,0,2", null},
                { "TTAN", "TGCA", "2,0,0,2", null},
                { "NTTA", "TGCA", "3,0,0,1", null},
                { "NGGA", "TGCA", "1,2,0,1", null},
                { "NTGGA", "TGCA", "2,2,0,1", null}
        };

        return testData;
    }

    @Test(dataProvider = "testData")
    public void testBase2Key(final String basesAsString, final String flowOrder,
                                    final String expectedKeyAsString, final String expectedClippedKeyAsString) {

        // int version
        final int[]         intKey = FlowBasedKeyCodec.base2key(basesAsString.getBytes(), flowOrder);
        Assert.assertNotNull(intKey);
        final String        intKeyAsString = StringUtils.join(intKey, ',');
        Assert.assertEquals(intKeyAsString, expectedKeyAsString);

        // byte version
        final byte[]         byteKey = FlowBasedKeyCodec.base2key(basesAsString.getBytes(), flowOrder, 1000);
        Assert.assertNotNull(byteKey);
        final String        byteKeyAsString = StringUtils.join(byteKey, ',');
        Assert.assertEquals(byteKeyAsString,
                expectedClippedKeyAsString != null ? expectedClippedKeyAsString : expectedKeyAsString);
    }

    @DataProvider(name="makeReadArrayTests")
    public Object[][] makeByteArrayTests(){
        List<Object[]> tests = new ArrayList<>();
        tests.add(new Object[]{new byte[]{'T','T','T','A','T','G','C'}, new byte[]{10,10,10,10,10,10,10}, "ACTG", (byte)0, new byte[]{0,0,10,10,10,10,10,10,10,10}});
        tests.add(new Object[]{new byte[]{'T','T','T','A','T','G','C'}, new byte[]{10,10,10,10,10,10,10}, "ACTG", (byte)10, new byte[]{10,10,10,10,10,10,10,10,10,10}});
        tests.add(new Object[]{new byte[]{'T','T','T','A','T','G','C'}, new byte[]{10,5,10,10,10,10,10}, "ACTG", (byte)0, new byte[]{0,0,5,5,10,10,10,10,10,10}});
        tests.add(new Object[]{new byte[]{'T','T','T','A','T','G','C'}, new byte[]{10,25,10,10,10,10,10}, "ACTG", (byte)0, new byte[]{0,0,10,10,10,10,10,10,10,10}});
        tests.add(new Object[]{new byte[]{'T','T','T','A','T','G','C'}, new byte[]{1,2,3,4,5,6,7}, "ACTG", (byte)0, new byte[]{0,0,1,1,4,4,5,6,6,7}});

        return tests.toArray(new Object[][]{});
    }

    @Test (dataProvider = "makeReadArrayTests")
    public void testBaseArray2KeySpace(final byte[] readBases, final byte[] qualArray, final String flowOrder, final byte defualtQual, final byte[] expected) {
        final int[] flowBases = FlowBasedKeyCodec.base2key(readBases, flowOrder);

        final byte[] result = FlowBasedKeyCodec.baseArray2KeySpace(readBases, flowBases.length, qualArray, defualtQual, flowOrder);

        Assert.assertEquals(flowBases.length, result.length, "Read bases in flow space and baseArray2KeySpace do not match in length");
        Assert.assertEquals(result, expected);
    }

}
