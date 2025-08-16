package org.broadinstitute.hellbender.utils.read;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class FlowBasedKeyCodecUnitTest extends GATKBaseTest {

    @DataProvider(name = "testData")
    public Object[][] getTestData() {

        final Object[][]        testData = {

                // trivial cases
                { "T", "TGCA", "1", "0",false },
                { "TT", "TGCA", "2", "0,0", false },
                { "TGCA", "TGCA", "1,1,1,1", "0,1,2,3",false },
                { "TA", "TGCA", "1,0,0,1", "0,3",false },
                { "TTAATG", "TGCA", "2,0,0,2,1,1", "0,0,3,3,4,5", false },

                // clipping
                { "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
                          , "TGCA", "130", String.join(",", Collections.nCopies(130, "0")), false },

                // N processing
                { "TNTA", "TGCA", "3,0,0,1", "0,0,0,3", false},
                { "TTNA", "TGCA", "3,0,0,1", "0,0,0,3", false},
                { "TTAN", "TGCA", "2,0,0,2", "0,0,3,3",false},
                { "TTAN", "TGCA", "2,0,0,2", "0,0,3,3", false},
                { "NTTA", "TGCA", "3,0,0,1", "0,0,0,3", false},
                { "NGGA", "TGCA", "1,2,0,1", "0,1,1,3", false},
                { "NTGGA", "TGCA", "2,2,0,1", "0,0,1,1,3", false},

                { "NT*GGA", "TGCA", "2,2,0,1", "0,0,1,1,3", true}
        };

        return testData;
    }

    @Test(dataProvider = "testData")
    public void testBase2Key(final String basesAsString, final String flowOrder,
                                    final String expectedKeyAsString, final String unused, final boolean expectException) {

        // int version
        try {
            final int[] intKey = FlowBasedKeyCodec.baseArrayToKey(basesAsString.getBytes(), flowOrder);
            Assert.assertNotNull(intKey);
            final String        intKeyAsString = StringUtils.join(intKey, ',');
            Assert.assertEquals(intKeyAsString, expectedKeyAsString);
        } catch (Exception e) {
            if ( !expectException )
                throw e;
        }
    }

    @Test(dataProvider = "testData")
    public void testKey2Flow(final String basesAsString, final String flowOrder,
                             final String key,
                             final String expectedBase2Flow,
                             boolean expectException) {

        // int version
        try {
            final int[] intKey = Arrays.stream(StringUtils.split(key, ','))
                    .mapToInt(Integer::parseInt).toArray();
            final int[] baseToFlow = FlowBasedKeyCodec.getBaseToFlow(intKey);
            final int[] expectedBase2FlowInt = Arrays.stream(StringUtils.split(expectedBase2Flow, ','))
                    .mapToInt(Integer::parseInt).toArray();
            Assert.assertEquals(baseToFlow.length, expectedBase2FlowInt.length);
            Assert.assertEquals(baseToFlow, expectedBase2FlowInt);
        } catch (Exception e) {
            if ( !expectException )
                throw e;
        }
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
        final int[] flowBases = FlowBasedKeyCodec.baseArrayToKey(readBases, flowOrder);

        final byte[] result = FlowBasedKeyCodec.baseArrayToKeySpace(readBases, flowBases.length, qualArray, defualtQual, flowOrder);

        Assert.assertEquals(flowBases.length, result.length, "Read bases in flow space and baseArrayToKeySpace do not match in length");
        Assert.assertEquals(result, expected);
    }

}
