package org.ultimagenomics.variant_calling;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.ultimagenomics.flow_based_read.read.FlowBasedRead;
import org.ultimagenomics.flow_based_read.utils.FlowBasedAlignmentArgumentCollection;

import java.nio.channels.SeekableByteChannel;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Function;

public class TrimmedReadsReader {

    private static final Logger logger = LogManager.getLogger(TrimmedReadsReader.class);

    private SamReader               samReader;
    private CountingReadFilter      readFilter;
    private Map<String, Integer>    readGroupMaxClass = new LinkedHashMap<>();
    private Map<String, String>     readGroupFlowOrder = new LinkedHashMap<>();
    private FlowBasedAlignmentArgumentCollection fbArgs = new FlowBasedAlignmentArgumentCollection();

    public TrimmedReadsReader(HaplotypeBasedVariantRecallerArgumentCollection vrArgs, Path referencePath) {
        Path samPath = IOUtils.getPath(vrArgs.READS_BAM_FILE);

        Function<SeekableByteChannel, SeekableByteChannel> cloudWrapper = BucketUtils.getPrefetchingWrapper(40);
        Function<SeekableByteChannel, SeekableByteChannel> cloudIndexWrapper = BucketUtils.getPrefetchingWrapper(40);

        samReader = SamReaderFactory.makeDefault().referenceSequence(referencePath).open(samPath, cloudWrapper, cloudIndexWrapper);
    }

    SAMSequenceDictionary getSamSequenceDictionary() {
        return samReader.getFileHeader().getSequenceDictionary();
    }

    public Collection<FlowBasedRead>  getReads(Locatable span, Locatable vcLoc) {

        List<FlowBasedRead>     reads = new LinkedList<>();
        SAMRecordIterator       iter = samReader.query(span.getContig(), span.getStart(), span.getEnd(), false);
        while ( iter.hasNext() ) {

            // establish record. ignore if variant context is not covered by this read?
            SAMRecord       record = iter.next();
            if ( !record.contains(vcLoc) )
                continue;

            // convert to gatk read
            String          readGroup = record.getReadGroup().getId();
            GATKRead        gatkRead = new SAMRecordToGATKReadAdapter(record);

            // filter out?
            if ( readFilter != null && !readFilter.test(gatkRead) )
                continue;

            // soft/hard clipped bases
            gatkRead = ReadClipper.hardClipSoftClippedBases(gatkRead);
            gatkRead = ReadClipper.hardClipToRegion(gatkRead, span.getStart(), span.getEnd());

            // convert to a flow based read
            int             maxClass = getMaxClass(readGroup);
            String          flowOrder = getFlowOrder(readGroup);
            FlowBasedRead   fbr = new FlowBasedRead(gatkRead, flowOrder, maxClass, fbArgs);
            fbr.apply_alignment();

            // clip to given span
            int read_start = fbr.getStart();
            int read_end = fbr.getEnd();
            int diff_left = span.getStart() - read_start;
            int diff_right = read_end - span.getEnd();
            fbr.apply_base_clipping(Math.max(0, diff_left), Math.max(diff_right, 0));

            // check if read is valid. it is possible that read was made invalid by apply_base_clipping
            // if so, ignore it (see FlowBasedRead.java:478 valid_key=false;
            if ( !fbr.is_valid() )
                continue;

            // add to output collection
            reads.add(fbr);
        }
        iter.close();

        return reads;
    }

    private synchronized int getMaxClass(String rg) {

        Integer     v = readGroupMaxClass.get(rg);

        if ( v == null ) {
            String mc_string = samReader.getFileHeader().getReadGroup(rg).getAttribute("mc");
            if ( mc_string == null ) {
                v = 12;
            } else {
                v = Integer.parseInt(mc_string);
            }
            readGroupMaxClass.put(rg, v);
            readGroupFlowOrder.put(rg, samReader.getFileHeader().getReadGroup(rg).getFlowOrder());
        }

        return v;
    }

    private String getFlowOrder(String rg) {
        return readGroupFlowOrder.get(rg);
    }

    public SAMFileHeader getHeader() {
        return samReader.getFileHeader();
    }

    public void setReadFilter(CountingReadFilter readFilter) {
        this.readFilter = readFilter;
    }

}
