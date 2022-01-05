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
import org.ultimagen.flowBasedRead.read.FlowBasedRead;
import org.ultimagen.flowBasedRead.utils.FlowBasedAlignmentArgumentCollection;

import java.nio.channels.SeekableByteChannel;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Function;

public class TrimmedReadsReader {

    private static final Logger logger = LogManager.getLogger(TrimmedReadsReader.class);

    private List<SamReader>         samReaders = new LinkedList<>();
    private CountingReadFilter      readFilter;
    private Map<String, Integer>    readGroupMaxClass = new LinkedHashMap<>();
    private Map<String, String>     readGroupFlowOrder = new LinkedHashMap<>();
    private FlowBasedAlignmentArgumentCollection fbArgs = new FlowBasedAlignmentArgumentCollection();

    public TrimmedReadsReader(List<Path> readsFiles, Path referencePath, int cloudPrefetchBuffer) {

        Function<SeekableByteChannel, SeekableByteChannel> cloudWrapper = BucketUtils.getPrefetchingWrapper(cloudPrefetchBuffer);
        Function<SeekableByteChannel, SeekableByteChannel> cloudIndexWrapper = BucketUtils.getPrefetchingWrapper(cloudPrefetchBuffer);

        for ( Path readsFile : readsFiles )
            samReaders.add(SamReaderFactory.makeDefault().referenceSequence(referencePath).open(readsFile, cloudWrapper, cloudIndexWrapper));
    }

    SAMSequenceDictionary getSamSequenceDictionary(SamReader samReader) {
        if ( samReader == null )
            samReader = samReaders.get(0);
        return samReader.getFileHeader().getSequenceDictionary();
    }

    public Map<SamReader, Collection<FlowBasedRead>>  getReads(Locatable span, Locatable vcLoc) {

        Map<SamReader, Collection<FlowBasedRead>>   readsByReader = new LinkedHashMap<>();
        for ( SamReader samReader : samReaders ) {
            List<FlowBasedRead>     reads = new LinkedList<>();
            SAMRecordIterator iter = samReader.query(span.getContig(), span.getStart(), span.getEnd(), false);
            while (iter.hasNext()) {

                // establish record. ignore if variant context is not covered by this read?
                SAMRecord record = iter.next();
                if (!record.contains(vcLoc))
                    continue;

                // convert to gatk read
                String readGroup = record.getReadGroup().getId();
                GATKRead gatkRead = new SAMRecordToGATKReadAdapter(record);

                // filter out?
                if (readFilter != null && !readFilter.test(gatkRead))
                    continue;

                // soft/hard clipped bases
                gatkRead = ReadClipper.hardClipSoftClippedBases(gatkRead);
                gatkRead = ReadClipper.hardClipToRegion(gatkRead, span.getStart(), span.getEnd());
                if (gatkRead.isUnmapped() || gatkRead.getCigar().isEmpty())
                    continue;

                // convert to a flow based read
                int maxClass = getMaxClass(readGroup, samReader.getFileHeader());
                String flowOrder = getFlowOrder(readGroup);
                FlowBasedRead fbr = new FlowBasedRead(gatkRead, flowOrder, maxClass, fbArgs);
                fbr.applyAlignment();

                // clip to given span
                int read_start = fbr.getStart();
                int read_end = fbr.getEnd();
                int diff_left = span.getStart() - read_start;
                int diff_right = read_end - span.getEnd();
                fbr.applyBaseClipping(Math.max(0, diff_left), Math.max(diff_right, 0), true);

                // check if read is valid. it is possible that read was made invalid by applyBaseClipping
                // if so, ignore it (see FlowBasedRead.java:478 valid_key=false;
                if (!fbr.isValid())
                    continue;

                // add to output collection
                reads.add(fbr);
            }
            iter.close();
            readsByReader.put(samReader, reads);
        }

        return readsByReader;
    }

    private synchronized int getMaxClass(String rg, SAMFileHeader hdr) {

        Integer     v = readGroupMaxClass.get(rg);

        if ( v == null ) {
            String mc_string = hdr.getReadGroup(rg).getAttribute("mc");
            if ( mc_string == null ) {
                v = 12;
            } else {
                v = Integer.parseInt(mc_string);
            }
            readGroupMaxClass.put(rg, v);
            readGroupFlowOrder.put(rg, hdr.getReadGroup(rg).getFlowOrder());
        }

        return v;
    }

    private String getFlowOrder(String rg) {
        return readGroupFlowOrder.get(rg);
    }

    public SAMFileHeader getHeader(SamReader samReader) {
        if ( samReader == null )
            samReader = samReaders.get(0);
        return samReader.getFileHeader();
    }

    public void setReadFilter(CountingReadFilter readFilter) {
        this.readFilter = readFilter;
    }

}
