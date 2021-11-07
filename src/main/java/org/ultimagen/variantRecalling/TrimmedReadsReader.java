package org.ultimagen.variantRecalling;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.ultimagen.FlowConstants;
import org.ultimagen.flowBasedRead.read.FlowBasedRead;
import org.ultimagen.flowBasedRead.utils.FlowBasedAlignmentArgumentCollection;

import java.nio.channels.SeekableByteChannel;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Function;

/**
 * a service class for HaplotypeBasedVariableRecaller that reads SAM/BAM files.
 *
 * For each given location (query location) it returns reads from all files that fall into the region after
 * trimming them accordign to span of a given variant context
 */
public class TrimmedReadsReader {

    private final List<SamReader>         samReaders = new LinkedList<>();
    private CountingReadFilter            readFilter;
    private final Map<String, Integer>    readGroupMaxClass = new LinkedHashMap<>();
    private final Map<String, String>     readGroupFlowOrder = new LinkedHashMap<>();
    private final FlowBasedAlignmentArgumentCollection fbArgs = new FlowBasedAlignmentArgumentCollection();

    public TrimmedReadsReader(final List<Path> readsFiles, final Path referencePath, final int cloudPrefetchBuffer) {

        final Function<SeekableByteChannel, SeekableByteChannel> cloudWrapper = BucketUtils.getPrefetchingWrapper(cloudPrefetchBuffer);
        final Function<SeekableByteChannel, SeekableByteChannel> cloudIndexWrapper = BucketUtils.getPrefetchingWrapper(cloudPrefetchBuffer);

        for ( Path readsFile : readsFiles ) {
            samReaders.add(SamReaderFactory.makeDefault().referenceSequence(referencePath).open(readsFile, cloudWrapper, cloudIndexWrapper));
        }
    }

    protected SAMSequenceDictionary getSamSequenceDictionary(final SamReader samReader) {
        return ((samReader != null) ? samReader : samReaders.get(0)).getFileHeader().getSequenceDictionary();
    }

    protected Map<SamReader, Collection<FlowBasedRead>>  getReads(final Locatable span, final Locatable vcLoc) {

        final Map<SamReader, Collection<FlowBasedRead>>   readsByReader = new LinkedHashMap<>();
        for ( SamReader samReader : samReaders ) {
            final List<FlowBasedRead>     reads = new LinkedList<>();
            final SAMRecordIterator iter = samReader.query(span.getContig(), span.getStart(), span.getEnd(), false);
            while (iter.hasNext()) {

                // establish record. ignore if variant context is not covered by this read?
                SAMRecord record = iter.next();
                if (!record.contains(vcLoc)) {
                    continue;
                }

                // convert to gatk read
                final String readGroup = record.getReadGroup().getId();
                GATKRead gatkRead = new SAMRecordToGATKReadAdapter(record);

                // filter out?
                if (readFilter != null && !readFilter.test(gatkRead)) {
                    continue;
                }

                // soft/hard clipped bases
                gatkRead = ReadClipper.hardClipSoftClippedBases(gatkRead);
                gatkRead = ReadClipper.hardClipToRegion(gatkRead, span.getStart(), span.getEnd());
                if (gatkRead.isUnmapped() || gatkRead.getCigar().isEmpty())
                    continue;

                // convert to a flow based read
                final int maxClass = getMaxClass(readGroup, samReader.getFileHeader());
                final String flowOrder = getFlowOrder(readGroup);
                final FlowBasedRead fbr = new FlowBasedRead(gatkRead, flowOrder, maxClass, fbArgs);
                fbr.applyAlignment();

                // clip to given span
                final int read_start = fbr.getStart();
                final int read_end = fbr.getEnd();
                final int diff_left = span.getStart() - read_start;
                final int diff_right = read_end - span.getEnd();
                fbr.applyBaseClipping(Math.max(0, diff_left), Math.max(diff_right, 0), true);

                // check if read is valid. it is possible that read was made invalid by applyBaseClipping
                // if so, ignore it (see FlowBasedRead.java:478 valid_key=false;
                if (!fbr.isValid()) {
                    continue;
                }

                // add to output collection
                reads.add(fbr);
            }
            iter.close();
            readsByReader.put(samReader, reads);
        }

        return readsByReader;
    }

    private synchronized int getMaxClass(final String rg, final SAMFileHeader hdr) {

        Integer     v = readGroupMaxClass.get(rg);

        if ( v == null ) {
            String mc = hdr.getReadGroup(rg).getAttribute("mc");
            if ( mc == null ) {
                v = FlowConstants.MAX_CLASS;
            } else {
                v = Integer.parseInt(mc);
            }
            readGroupMaxClass.put(rg, v);
            readGroupFlowOrder.put(rg, hdr.getReadGroup(rg).getFlowOrder());
        }

        return v;
    }

    private String getFlowOrder(final String rg) {
        return readGroupFlowOrder.get(rg);
    }

    public SAMFileHeader getHeader(final SamReader samReader) {
        return ((samReader != null) ? samReader : samReaders.get(0)).getFileHeader();
    }

    public void setReadFilter(final CountingReadFilter readFilter) {
        this.readFilter = readFilter;
    }

}
