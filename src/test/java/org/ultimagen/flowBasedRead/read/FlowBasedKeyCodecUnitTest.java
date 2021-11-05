package org.ultimagen.flowBasedRead.read;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.ultimagen.variantRecalling.TrimmedReadsReaderUnitTest;

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
}
