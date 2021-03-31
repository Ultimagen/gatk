package org.ultimagen.variantRecalling;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

public class VariantRecallerResultWriterTest extends GATKBaseTest {

    @Test
    public void testGetOffsetOnRead() throws Exception {

        // create sam header
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader(1, 1, 1000);

        // create reads
        final GATKRead read2 = ArtificialReadUtils.createArtificialRead("10M1D10M1I10M");
        int[]          ofsPairs = { 0, 0, 9, 9,
                                    10, 9,
                                    11, 10, 20, 19,
                                    21, 20,
                                    22, 22, 30, 30,
                                    31, -1
                                    };
        for ( int n = 0 ; n < ofsPairs.length ; n += 2 )
            Assert.assertEquals(VariantRecallerResultWriter.getOffsetOnRead(read2, ofsPairs[n]), ofsPairs[n+1]);


    }
}
