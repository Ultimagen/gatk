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
import picard.sam.markduplicates.util.ReadEnds;

import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Class representing a single read fragment at a particular start location without a mapped mate.
 *
 * This class holds onto as little information as possible in an attempt to prevent excessive serialization of
 * during the processing step of MarkDuplicatesSpark
 */
public class FlowFragment extends Fragment {

    protected final short flowScore;

    private final static Logger logger = LogManager.getLogger(FlowFragment.class);

    public FlowFragment(final GATKRead first, final SAMFileHeader header, int partitionIndex, MarkDuplicatesScoringStrategy scoringStrategy, Map<String, Byte> headerLibraryMap, final MarkDuplicatesSparkArgumentCollection mdArgs) {
        super(first, header, partitionIndex, scoringStrategy, headerLibraryMap, mdArgs);

        int        start = first.isReverseStrand() ? ReadUtils.getSelectedRecordEnd(first, null, header, mdArgs) : ReadUtils.getSelectedRecordStart(first, null, header, mdArgs);
        int        endUncert = 0;
        if ( mdArgs.FLOW_END_LOCATION_SIGNIFICANT ) {
            AtomicInteger endUncertainty = new AtomicInteger(mdArgs.ENDS_READ_UNCERTAINTY);
            this.end = !first.isReverseStrand() ? ReadUtils.getSelectedRecordEnd(first, endUncertainty, header, mdArgs) : ReadUtils.getSelectedRecordStart(first, endUncertainty, header, mdArgs);
            endUncert = endUncertainty.intValue();
        }
        this.key = ReadsKey.getKeyForFragment(start,
                isRead1ReverseStrand(),
                (short)ReadUtils.getReferenceIndex(first, header),
                headerLibraryMap.get(MarkDuplicatesSparkUtils.getLibraryForRead(first, header, LibraryIdGenerator.UNKNOWN_LIBRARY)),
                mdArgs);

        this.flowScore = (this.end != ReadUtils.FLOW_BASED_INSIGNIFICANT_END_UNCERTIANTY)
                ? ((mdArgs.FLOW_QUALITY_SUM_STRATEGY && isFlow(first)) ? computeFlowDuplicateScore(first, start, end) : scoringStrategy.score(first))
                : -1;
        if ( mdArgs.DEBUG_ULTIMA_DUPS || isDebugUltimaRead(mdArgs, first) ) {
            logger.info(String.format("F [%s %s] : unc %d-%d : clp %d-%d -> %d-%d(%d) %d",
                    first.getName(), first.isReverseStrand() ? "R" : "N",
                    first.getUnclippedStart(), first.getUnclippedEnd(),
                    first.getStart(), first.getEnd(),
                    start, end, endUncert, this.flowScore));
        }

    }

    private boolean isDebugUltimaRead(MarkDuplicatesSparkArgumentCollection mdArgs, GATKRead first) {
        return mdArgs.DEBUG_ULTIMA_READ_NAME != null && mdArgs.DEBUG_ULTIMA_READ_NAME.contains(first.getName());
    }

    private short computeFlowDuplicateScore(GATKRead rec, int start, int end) {

        Short storedScore = (Short)rec.getTransientAttribute("DuplicateScore");
        if ( storedScore == null ) {
            short score = 0;

            score += (short) Math.min(getFlowSumOfBaseQualities(rec, start, end), Short.MAX_VALUE / 2);

            storedScore = score;
            rec.setTransientAttribute("DuplicateScore", storedScore);
        }

        return storedScore;
    }

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
