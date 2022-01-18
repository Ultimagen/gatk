package org.broadinstitute.hellbender.utils.read.markduplicates.sparkrecords;

import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.transforms.markduplicates.MarkDuplicatesSparkUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.markduplicates.LibraryIdGenerator;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import org.broadinstitute.hellbender.utils.read.markduplicates.ReadsKey;

import java.util.Map;

/**
 * Flow-based specific class  representing a single read fragment at a particular start location without a mapped mate.
 *
 * This class holds onto as little information as possible in an attempt to prevent excessive serialization of
 * during the processing step of MarkDuplicatesSpark
 *
 * It differs from its parent class in that it takes into account the end of the fragment as well (if instructed
 * to by FLOW_END_LOCATION_SIGNIFICANT). It also allows for an alternative flow-specific score function to be deployed.
 */
public class FlowFragment extends Fragment {
    static final long serialVersionUID = 1L;
    public static final String FLOW_DUPLICATE_SCORE_ATTR_NAME = "FlowDuplicateScore";
    protected final short flowScore;

    public FlowFragment(final GATKRead first, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy, Map<String, Byte> headerLibraryMap, final MarkDuplicatesSparkArgumentCollection mdArgs) {
        super(first, header, partitionIndex, scoringStrategy, headerLibraryMap, mdArgs);

        int        start = first.isReverseStrand() ? ReadUtils.getMarkDupReadEnd(first, false, header, mdArgs) : ReadUtils.getMarkDupReadStart(first, false, header, mdArgs);
        if ( mdArgs.FLOW_END_LOCATION_SIGNIFICANT ) {
            this.end = !first.isReverseStrand() ? ReadUtils.getMarkDupReadEnd(first, true, header, mdArgs) : ReadUtils.getMarkDupReadStart(first, true, header, mdArgs);
        }
        this.key = ReadsKey.getKeyForFragment(start,
                isRead1ReverseStrand(),
                (short)ReadUtils.getReferenceIndex(first, header),
                headerLibraryMap.get(MarkDuplicatesSparkUtils.getLibraryForRead(first, header, LibraryIdGenerator.UNKNOWN_LIBRARY)),
                mdArgs);

        this.flowScore = (this.end != ReadUtils.FLOW_BASED_INSIGNIFICANT_END)
                ? ((mdArgs.FLOW_QUALITY_SUM_STRATEGY && isFlow(first)) ? computeFlowDuplicateScore(first, start, end) : scoringStrategy.score(first))
                : -1;
    }

    // compute fragment score using a flow-based specific method - cache in transient attribute
    private short computeFlowDuplicateScore(GATKRead rec, int start, int end) {

        Short storedScore = (Short)rec.getTransientAttribute(FLOW_DUPLICATE_SCORE_ATTR_NAME);
        if ( storedScore == null ) {
            short score = 0;

            score += (short) Math.min(getFlowSumOfBaseQualities(rec, start, end), Short.MAX_VALUE / 2);

            storedScore = score;
            rec.setTransientAttribute(FLOW_DUPLICATE_SCORE_ATTR_NAME, storedScore);
        }

        return storedScore;
    }

    /**
     * A quality summing scoring strategy used for flow based reads.
     *
     * The method walks on the bases of the read, in the synthesis direction. For each base, the effective
     * quality value is defined as the value on the first base on the hmer to which the base belongs to. The score
     * is defined to be the sum of all effective values above a given threshold.
     */
    private int getFlowSumOfBaseQualities(GATKRead rec, int start, int end) {
        int score = 0;

        if ( rec.isReverseStrand() ) {
            int     tmp = start;
            start = end;
            end = tmp;
        }

        // access qualities and bases
        byte[]      quals = rec.getBaseQualitiesNoCopy();
        byte[]      bases = rec.getBasesNoCopy();

        // establish range of bases/quals to work on
        int         startingOffset = Math.max(0, start - rec.getUnclippedStart());
        int         endOffset = Math.max(0, rec.getUnclippedEnd() - end);

        // loop on bases, extract qual related to homopolymer from start of homopolymer
        byte        lastBase = 0;
        byte        effectiveQual = 0;
        for ( int i = startingOffset ; i < bases.length - endOffset ; i++ ) {
            byte        base = bases[i];
            if ( base != lastBase ) {
                effectiveQual = quals[i];
            }
            if ( effectiveQual >= 15 ) {
                score += effectiveQual;
            }
            lastBase = base;
        }

        return score;
    }

    private static boolean isFlow(GATKRead rec) {
        return rec.hasAttribute("tp");
    }

    @Override
    public short getScore() {
        return flowScore;
    }

}
