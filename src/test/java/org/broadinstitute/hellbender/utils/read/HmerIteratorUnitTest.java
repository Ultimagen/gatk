package org.broadinstitute.hellbender.utils.read;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class HmerIteratorUnitTest extends GATKBaseTest {

    @DataProvider(name="DataProvider")
    public Object[][] getDataProvider() {
        return new Object[][] {
            {
                "AAAABBBCCD", new int[]{'A', 4, 'B', 3, 'C', 2, 'D', 1}
            }
        };
    }

    @Test(dataProvider = "DataProvider")
    public void testHmerIterator(final String bases, int[] expectedResults) {

        final BaseUtils.HmerIterator iter = new BaseUtils.HmerIterator(bases.getBytes());
        int i = 0;

        // loop over hmers
        while ( iter.hasNext() ) {
            Pair<Byte,Integer>      pair = iter.next();
            Assert.assertEquals(pair.getLeft().intValue(), expectedResults[i]);
            Assert.assertEquals(pair.getRight().intValue(), expectedResults[i+1]);
            i += 2;
        }

        // must be at end
        Assert.assertEquals(i, expectedResults.length);
    }
}
