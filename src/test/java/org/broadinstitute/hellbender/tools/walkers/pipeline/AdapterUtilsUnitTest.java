package org.broadinstitute.hellbender.tools.walkers.pipeline;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class AdapterUtilsUnitTest extends GATKBaseTest {

    @DataProvider(name = "iupacMatch")
    Object[][] iupacMatchDataProvider() {
        return new Object[][] {
                // iupa, seq, result
                { "A", "ACGTU", "10000" },
                { "C", "ACGTU", "01000" },
                { "G", "ACGTU", "00100" },
                { "T", "ACGTU", "00010" },
                { "U", "ACGTU", "00001" },

                { "X", "ACGTU", "11110" },
                { "N", "ACGTU", "11110" },

                { "M", "ACGT", "1100" },
                { "R", "ACGT", "1010" },
                { "W", "ACGT", "1001" },
                { "S", "ACGT", "0110" },
                { "Y", "ACGT", "0101" },
                { "K", "ACGT", "0011" },

                { "V", "ACGT", "1110" },
                { "H", "ACGT", "1101" },
                { "D", "ACGT", "1011" },
                { "B", "ACGT", "0111" },

                { "?", "ACGT", "0000" },
        };
    }

    @Test(dataProvider = "iupacMatch")
    public void testIupacMatch(final String pattern, final String seq, final String result) {
        final byte iupac = pattern.getBytes()[0];
        final byte[] bases = seq.getBytes();
        final StringBuilder sb = new StringBuilder();
        for ( int i = 0 ; i < bases.length ; i++ ) {
            sb.append(AdapterUtils.iupacMatch(bases[i], iupac) ? '1' : '0');
        }
        Assert.assertEquals(result, sb.toString());
    }

    @DataProvider(name = "findAdapter")
    public Object[][] findAdapterDataProvider() {
        return new Object[][] {
                // read, pattern, error_rate, min_overlap, found_start, found_length

                // trivial
                { "ACGT", "A", 0.0f, 10, 0, 1 },
                { "ACGT", "C", 0.0f, 10, 1, 1 },
                { "ACGT", "G", 0.0f, 10, 2, 1 },
                { "ACGT", "T", 0.0f, 10, 3, 1 },
                { "ACGT", "U", 0.0f, 10, -1, 0 },

                // multi byte
                { "ACGT", "AC", 0.0f, 10, 0, 2 },
                { "ACGT", "CGT", 0.0f, 10, 1, 3 },
                { "ACGT", "GT", 0.0f, 10, 2, 2 },

                // error rate
                { "ATTTTTTTTTTG", "TTTTTTTTTT", 0.0f, 0, 1, 10 },  // exact
                { "ATTTTTTTTTTG", "CTTTTTTTTT", 0.0f, 0, -1, 0 },  // exact fail
                { "ATTTTTTTTTTG", "CTTTTTTTTT", 0.1f, 0, 0, 10 },  // 0.1 error rate
                { "ATTTTTTTTTTG", "CCTTTTTTTT", 0.1f, 0, -1, 0 },  // 0.1 error rate fail
                { "ATTTTTTTTTTG", "TTTTTTTTTT", 0.1f, 0, 1, 10 },  // 0.1 error rate find best

                // min overlap (note that errorRate is set to 1.0 to allow any errors)
                { "ATTTTTTTTTTG", "TTTTTTTTTT", 1.0f, 10, 1, 10 },  // exact
                { "ATTTTTTTTTTG", "CTTTTTTTTT", 1.0f, 10, -1, 0 },  // exact, fail

                // anchors
                { "ACGT", "^", 0.0f, 10, 0, 0 },
                { "ACGT", "$", 0.0f, 10, 4, 0 },
                { "ATTTTTTTTTTG", "^TTTTTTTTTT", 0.0f, 0, -1, 0 },
                { "ATTTTTTTTTTG", "^ATTTTTTTTT", 0.0f, 0, 0, 10 },
                { "ATTTTTTTTTTG", "TTTTTTTTTT$", 0.0f, 0, -1, 0 },
                { "ATTTTTTTTTTG", "TTTTTTTTTG$", 0.0f, 0, 2, 10 },

        };
    }

    @Test(dataProvider = "findAdapter")
    public void testFindAdapter(String read, String pattern, double errorRate, int minOverlap, int foundStart, int foundLength) {

        AdapterUtils.AdapterPattern ap = new AdapterUtils.AdapterPattern(pattern, errorRate, minOverlap, false, false);
        AdapterUtils.FoundAdapter fa = AdapterUtils.findAdapter(read.getBytes(), ap, 0, read.getBytes().length);

        if ( fa == null ) {
            Assert.assertEquals(-1, foundStart, "wrong start. should have found " + pattern + " in " + read + " " + errorRate + "/" + minOverlap);
        } else {
            Assert.assertEquals(fa.start, foundStart, "wrong start on finding " + pattern + " in " + read + " " + errorRate + "/" + minOverlap);
            Assert.assertEquals(fa.length, foundLength, "wrong length on finding " + pattern + " in " + read + " " + errorRate + "/" + minOverlap);
        }
    }
}
