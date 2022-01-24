package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Lazy;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContextUtils;
import org.apache.commons.jexl2.Expression;
import org.broadinstitute.hellbender.tools.FlowBasedAlignmentArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.groundtruth.GroundTruthReadsBuilder;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * utility class for flow based read
 */
public class FlowBasedReadUtils {

    public static final int FLOW_SUM_OF_BASE_QUALITY_THRESHOLD = 15;

    private static final Map<String, ReadGroupInfo> readGroupInfo = new LinkedHashMap<>();

    static public class ReadGroupInfo {
        final public String  flowOrder;
        final public int     maxClass;

        private String  reversedFlowOrder = null;

        ReadGroupInfo(final SAMReadGroupRecord readGroup) {

            Utils.nonNull(readGroup);
            this.flowOrder = readGroup.getFlowOrder();
            Utils.nonNull(this.flowOrder);

            String mc = readGroup.getAttribute(FlowBasedRead.MAX_CLASS_READ_GROUP_TAG);
            this.maxClass = (mc == null) ? FlowBasedRead.MAX_CLASS : Integer.parseInt(mc);
        }

        public synchronized String getReversedFlowOrder() {
            if ( reversedFlowOrder == null ) {
                reversedFlowOrder = SequenceUtil.reverseComplement(flowOrder);
            }
            return reversedFlowOrder;
        }
    }

    public static boolean readEndMarkedUncertain(final GATKRead rec) {
        final String        tm = rec.getAttributeAsString(FlowBasedRead.CLIPPING_TAG_NAME);
        if ( tm == null ) {
            return false;
        } else {
            return tm.indexOf('Q') >= 0 || tm.indexOf('Z') >= 0;
        }
    }

    public static boolean readEndMarkedUnclipped(final GATKRead rec, boolean FLOW_Q_IS_KNOWN_END) {
        final String        tm = rec.getAttributeAsString(FlowBasedRead.CLIPPING_TAG_NAME);
        if ( tm == null ) {
            return false;
        } else {
            return (tm.indexOf('A') >= 0) || (FLOW_Q_IS_KNOWN_END && (tm.indexOf('Q') >= 0));
        }
    }

    // get flow order for a specific read
    public static byte[] getReadFlowOrder(final SAMFileHeader header, GATKRead read) {

        // are we looking for a specific read group, as specified by the read?
        final String    readGroupName = (read != null) ? read.getReadGroup() : null;
        if ( readGroupName != null ) {
            final SAMReadGroupRecord rg = header.getReadGroup(readGroupName);
            if ( rg != null && rg.getFlowOrder() != null )
                return rg.getFlowOrder().getBytes();
        }

        // if here, either no read was specified, or the read has no group, or the group is not found, or it has no flow
        // revert to old behavior of returning the first found
        for ( SAMReadGroupRecord rg : header.getReadGroups() ) {
            // must match read group name?
            String      flowOrder = rg.getFlowOrder();
            if ( flowOrder != null ) {
                return flowOrder.getBytes();
            }
        }
        return null;
    }

    /**
     * Computes the sum of base qualities of the given flow read.
     */
    public static int flowSumOfBaseQualities(final GATKRead read) {
        if (read == null) {
            return 0;
        } else {
            int sum = 0;

            // access qualities and bases
            byte[]      quals = read.getBaseQualitiesNoCopy();
            byte[]      bases = read.getBasesNoCopy();

            // loop on bases, extract qual related to homopolymer from start of homopolymer
            int         i = 0;
            byte        lastBase = 0;
            byte        effectiveQual = 0;
            for (final byte base : bases ) {
                if ( base != lastBase )
                    effectiveQual = quals[i];
                if ( effectiveQual >= FLOW_SUM_OF_BASE_QUALITY_THRESHOLD )
                    sum += effectiveQual;
                lastBase = base;
                i++;
            }

            return sum;
        }
    }

    public static boolean isFlow(final GATKRead rec) {
        return rec.hasAttribute(FlowBasedRead.FLOW_MATRIX_TAG_NAME)
                || rec.hasAttribute(FlowBasedRead.FLOW_MATRiX_OLD_TAG_KR)
                || rec.hasAttribute(FlowBasedRead.FLOW_MATRiX_OLD_TAG_TI);
    }

    public static boolean isFlow(final SAMRecord rec) {
        return rec.hasAttribute(FlowBasedRead.FLOW_MATRIX_TAG_NAME)
                || rec.hasAttribute(FlowBasedRead.FLOW_MATRiX_OLD_TAG_KR)
                || rec.hasAttribute(FlowBasedRead.FLOW_MATRiX_OLD_TAG_TI);

    }

    public static synchronized ReadGroupInfo getReadGroupInfo(final SAMFileHeader hdr, final GATKRead read) {

        String              name = read.getReadGroup();
        Utils.nonNull(name);
        ReadGroupInfo info = readGroupInfo.get(name);
        if ( info == null ) {
            readGroupInfo.put(name, info = new ReadGroupInfo(hdr.getReadGroup(name)));
        }
        return info;
    }

    /**
     * Finds  a usable FlowOrder to be used for engine calculation (when no specufic flow order already established for a specific read)
     */
    public static String findFirstUsableFlowOrder(final SAMFileHeader hdr, final FlowBasedAlignmentArgumentCollection fbargs) {
        for ( final SAMReadGroupRecord rg : hdr.getReadGroups() ) {
            final String flowOrder = rg.getFlowOrder();
            if ( flowOrder != null && flowOrder.length() >= fbargs.flowOrderCycleLength ) {
                return flowOrder.substring(0, fbargs.flowOrderCycleLength);
            }
        }
        return null;
    }

}
