package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * utility class for flow based read
 */
public class FlowBasedReadUtils {

    public static final int FLOW_SUM_OF_BASE_QUALITY_THRESHOLD = 15;
    public static final FlowBasedArgumentCollection DEFAULT_FLOW_BASED_ARGUMENT_COLLECTION = new FlowBasedArgumentCollection();

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

        if ( !isFlow(read) ) {
            throw new IllegalArgumentException("read must be flow based: " + read);
        }

        String              name = read.getReadGroup();
        Utils.nonNull(name);
        ReadGroupInfo info = readGroupInfo.get(name);
        if ( info == null ) {
            readGroupInfo.put(name, info = new ReadGroupInfo(hdr.getReadGroup(name)));
        }
        return info;
    }

    /**
     * Finds a usable FlowOrder to be used for engine calculation (when no specufic flow order already established for a specific read)
     */
    public static String findFirstUsableFlowOrder(final SAMFileHeader hdr, final FlowBasedArgumentCollection fbargs) {
        for ( final SAMReadGroupRecord rg : hdr.getReadGroups() ) {
            final String flowOrder = rg.getFlowOrder();
            if ( flowOrder != null && flowOrder.length() >= fbargs.flowOrderCycleLength ) {
                return flowOrder.substring(0, fbargs.flowOrderCycleLength);
            }
        }

        throw new GATKException("Unable to perform flow based operations without the flow order");
    }

    /*
     * clips flows from the left to clip the input number of bases
     * Needed to trim the haplotype to the read
     * Returns number of flows to remove and the change in the left most remaining flow if necessary
     */
    static public int[] findLeftClipping(final int baseClipping, final int[] flow2base, final int[] key) {
        final int [] result = new int[2];
        if (baseClipping == 0 ){
            return result;
        }

        int stopClip = 0;
        for (int i = 0 ; i < flow2base.length; i++ ) {

            if (flow2base[i] + key[i] >= baseClipping) {
                stopClip = i;
                break;
            }
        }
        final int hmerClipped = baseClipping - flow2base[stopClip] - 1;
        result[0] = stopClip;
        result[1] = hmerClipped;
        return result;
    }

    /*
     * clips flows from the right to trim the input number of bases
     * Returns number of flows to remove and the change in the right most flow.
     */
    static public int[] findRightClipping(final int baseClipping, final int[] rFlow2Base, final int[] rKey) {
        final int [] result = new int[2];
        if (baseClipping == 0 ){
            return result;
        }

        int stopClip = 0;

        for (int i = 0; i < rFlow2Base.length; i++ ) {
            if (rFlow2Base[i] + rKey[i] >= baseClipping) {
                stopClip = i;
                break;
            }
        }

        final int hmerClipped = baseClipping - rFlow2Base[stopClip] - 1;
        result[0] = stopClip;
        result[1] = hmerClipped;
        return result;
    }

    /**
     * create a FlowBasedRead from a proper SAMRecord
     */
    static public FlowBasedRead convertToFlowBasedRead(GATKRead read, SAMFileHeader header) {
        final ReadGroupInfo readGroupInfo = getReadGroupInfo(header, read);
        return new FlowBasedRead(read, readGroupInfo.flowOrder, readGroupInfo.maxClass, DEFAULT_FLOW_BASED_ARGUMENT_COLLECTION);
    }
}
