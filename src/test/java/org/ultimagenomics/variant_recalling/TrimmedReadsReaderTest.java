package org.ultimagenomics.variant_recalling;

import htsjdk.samtools.SamReader;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.ultimagenomics.flow_based_read.read.FlowBasedRead;

import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class TrimmedReadsReaderTest extends GATKBaseTest {

    protected static String    vcTestDir = publicTestDir + Const.DATA_DIR;

    @Test
    public void testBasic() throws Exception {

        // establish paths
        Path        path1 = FileSystems.getDefault().getPath(vcTestDir, "chr5.bam1.rename.bam");
        Path        path2 = FileSystems.getDefault().getPath(vcTestDir, "chr5.bam2.rename.bam");
        List<Path>  paths = new LinkedList<>();
        paths.add(path1);
        paths.add(path2);

        // create reader
        TrimmedReadsReader  reader = new TrimmedReadsReader(paths, null, 40);

        // global
        Assert.assertNotNull(reader.getSamSequenceDictionary(null));
        Assert.assertNotNull(reader.getHeader(null));

        // reads
        Map<SamReader, Collection<FlowBasedRead>> reads = reader.getReads(
                                new SimpleInterval("chr5:70036483-70036764"),
                                new SimpleInterval("chr5:70036625-70036625"));
        Assert.assertEquals(reads.size(), 2);
        for ( Map.Entry<SamReader, Collection<FlowBasedRead>> entry : reads.entrySet() ) {

            Assert.assertNotNull(reader.getSamSequenceDictionary(entry.getKey()));
            Assert.assertNotNull(reader.getHeader(entry.getKey()));

            Assert.assertNotEquals(entry.getValue().size(), 0);
        }
    }
}
