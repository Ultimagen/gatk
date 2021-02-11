package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.LibraryIdGenerator;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;
import picard.sam.markduplicates.util.ReadEnds;

import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Dummy class representing a mated read fragment at a particular start position to be used for accounting
 * when deciding whether to duplicate unmatched fragments.
 *
 * This class holds onto as little information as possible in an attempt to prevent excessive serialization of
 * during the processing step of MarkDuplicatesSpark
 */
public final class EmptyFragment extends PairedEnds {

    private final static Logger logger = LogManager.getLogger(EmptyFragment.class);

    protected transient ReadsKey key;

    private final boolean R1R;

    private final short score;

    /**
     * special constructor for creating empty fragments
     * this only includes the necessary data to locate the read, the rest is unnecessary because it will appear in the paired bucket
     *
     */
    public EmptyFragment(GATKRead read, SAMFileHeader header, Map<String, Byte> headerLibraryMap, final MarkDuplicatesSparkArgumentCollection mdArgs) {
        super(0, null);

        int        start = read.isReverseStrand() ? ReadUtils.getSelectedRecordEnd(read, null, header, mdArgs) : ReadUtils.getSelectedRecordStart(read, null, header, mdArgs);
        this.R1R = read.isReverseStrand();
        this.key = ReadsKey.getKeyForFragment(start,
                isRead1ReverseStrand(),
                ReadUtils.getReferenceIndex(read, header),
                headerLibraryMap.get(MarkDuplicatesSparkUtils.getLibraryForRead(read, header, LibraryIdGenerator.UNKNOWN_LIBRARY)),
                mdArgs);

        score = 0;
    }

    @Override
    public Type getType() {
        return Type.EMPTY_FRAGMENT;
    }
    @Override
    public short getScore() {
        return score;
    }
    @Override
    // NOTE: This is transient and thus may not exist if the object gets serialized
    public ReadsKey key() {
        return key;
    }
    @Override
    public boolean isRead1ReverseStrand() {
        return R1R;
    }
    @Override
    public byte getOrientationForPCRDuplicates() {
        return (R1R)? ReadEnds.R : ReadEnds.F;
    }
    @Override
    public String toString() {
        return "EmptyFragment ";
    }
}
