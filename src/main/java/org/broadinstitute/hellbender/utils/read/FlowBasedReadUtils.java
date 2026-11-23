package org.broadinstitute.hellbender.utils.read;

import gov.nih.nlm.ncbi.ngs.NGS;
import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.utils.NGSPlatform;
import org.broadinstitute.hellbender.utils.SequencerFlowClass;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;


/**
 * Utility class for working with flow-based reads
 *
 * The main member static class is {@code ReadGroupInfo} that contains methods that allow
 * working with headers of flow based reads and extracting the flow order and the maximal hmer class called.
 * It also contains methods to check how the read was clipped ({@code readEndMarkedUnclipped}, {@code readEndMarkedUncertain})
 * Lastly, {@code FlowBasedReadUtils.isFlowPlatform} is **the** function to determine if the data are flow-based
 *
 */
public class FlowBasedReadUtils {

    public static final int FLOW_SUM_OF_BASE_QUALITY_THRESHOLD = 15;
    public static final FlowBasedArgumentCollection DEFAULT_FLOW_BASED_ARGUMENT_COLLECTION = new FlowBasedArgumentCollection();
    static final public int FLOW_BASED_INSIGNIFICANT_END = 0;

    private static final Map<String, ReadGroupInfo> readGroupInfo = new LinkedHashMap<>();

    public enum CycleSkipStatus {
        NS(0),         // no skip
        PCS(1),        // possible cycle skip
        CS(2);         // cycle skip

        private int priority;

        CycleSkipStatus(final int priority) {
            this.priority = priority;
        }

        int getPriority() {
            return this.priority;
        }

    }

    static public class ReadGroupInfo {
        final public String  flowOrder;
        final public int     maxClass;
        final public boolean isFlowPlatform;
        private String  reversedFlowOrder = null;

        public ReadGroupInfo(final SAMReadGroupRecord readGroup) {
            if (readGroup.getPlatform()==null){
                isFlowPlatform = false;
            } else if (NGSPlatform.fromReadGroupPL(readGroup.getPlatform())==NGSPlatform.UNKNOWN){
                isFlowPlatform = false;
            } else if (NGSPlatform.fromReadGroupPL(readGroup.getPlatform()) == NGSPlatform.LS454) {
                //old Ultima data can have PLATFORM==LS454 and not have FO tag
                isFlowPlatform = true;
            } else if (NGSPlatform.fromReadGroupPL(readGroup.getPlatform()) == NGSPlatform.ULTIMA){
                if (readGroup.getFlowOrder()!=null) {
                    isFlowPlatform = true;
                } else {
                    throw new RuntimeException("Malformed Ultima read group identified, aborting: " + readGroup);
                }
            } else {
                isFlowPlatform = false;
            }

            if (isFlowPlatform) {
                this.flowOrder = readGroup.getFlowOrder();
                String mc = readGroup.getAttribute(FlowBasedRead.MAX_CLASS_READ_GROUP_TAG);
                this.maxClass = (mc == null) ? FlowBasedRead.MAX_CLASS : Integer.parseInt(mc);
            } else { // not a flow platform
                this.flowOrder = null;
                this.maxClass = 0;
            }
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

    /** name of the read group attribute holding the maximal QUAL character (tp:Z sparse format) */
    public static final String READ_GROUP_MAX_QUAL_ATTRIBUTE = "mq";

    /**
     * @return the maximal quality (phred, not ascii) declared by the read group's mq attribute, or null if absent
     */
    public static Byte getReadGroupMaxQual(final SAMReadGroupRecord rg) {
        if ( rg == null ) {
            return null;
        }
        final String mq = rg.getAttribute(READ_GROUP_MAX_QUAL_ATTRIBUTE);
        if ( mq == null ) {
            return null;
        }
        if ( mq.length() != 1 ) {
            throw new GATKException("Illegal read group mq attribute (expected a single character): " + mq);
        }
        return (byte)(mq.charAt(0) - 33);
    }

    /**
     * Decode the tp tag of a record into a dense per-base array of -1/0/+1.
     * Supports the old format (tp:B:c, one value per base) and the sparse format (tp:Z containing only '+' and '-'
     * for bases whose QUAL differs from the read group's mq; an absent tag means all zeros when mq is present).
     *
     * @param rec record
     * @param header header used to find the read group (if null, the record's own header is used)
     * @return dense tp array, or null if the read has no tp tag and its read group has no mq
     * @throws GATKException if the sparse tag is malformed, or mq is missing for a tp:Z tag
     */
    public static byte[] getTpValues(final SAMRecord rec, final SAMFileHeader header) {
        final Object tp = rec.getAttribute(FlowBasedRead.FLOW_MATRIX_TAG_NAME);
        if ( tp != null && !(tp instanceof String) ) {
            return rec.getSignedByteArrayAttribute(FlowBasedRead.FLOW_MATRIX_TAG_NAME);
        }

        final SAMFileHeader hdr = (header != null) ? header : rec.getHeader();
        final String rgName = (String)rec.getAttribute("RG");
        final Byte mq = (hdr != null && rgName != null) ? getReadGroupMaxQual(hdr.getReadGroup(rgName)) : null;
        if ( mq == null ) {
            if ( tp == null ) {
                return null;
            }
            throw new GATKException("tp:Z tag requires read group mq attribute: " + rec.getReadName());
        }

        final byte[] quals = rec.getBaseQualities();
        final byte[] result = new byte[quals.length];
        final String sparse = (tp == null) ? "" : (String)tp;
        int next = 0;
        for ( int i = 0 ; i < quals.length ; i++ ) {
            if ( quals[i] == mq ) {
                continue;
            }
            if ( next >= sparse.length() ) {
                throw new GATKException("tp string too short for read " + rec.getReadName());
            }
            final char c = sparse.charAt(next++);
            if ( c == '+' ) {
                result[i] = 1;
            } else if ( c == '-' ) {
                result[i] = -1;
            } else {
                throw new GATKException("Illegal character '" + c + "' in tp string of read " + rec.getReadName());
            }
        }
        if ( next != sparse.length() ) {
            throw new GATKException("tp string has leftover characters for read " + rec.getReadName());
        }
        return result;
    }

    /** same as {@link #getTpValues(SAMRecord, SAMFileHeader)} for a GATKRead */
    public static byte[] getTpValues(final GATKRead read, final SAMFileHeader header) {
        if ( read instanceof SAMRecordToGATKReadAdapter ) {
            return getTpValues(((SAMRecordToGATKReadAdapter)read).getEncapsulatedSamRecord(), header);
        }
        // no access to the underlying record type: only the old format can be handled
        return read.getAttributeAsByteArray(FlowBasedRead.FLOW_MATRIX_TAG_NAME);
    }

    /** true if the read has flow tags, or (header permitting) belongs to a read group declaring mq (absent tp == all zeros) */
    public static boolean hasFlowTags(final GATKRead rec, final SAMFileHeader hdr) {
        return hasFlowTags(rec) || (hdr != null && rec.getReadGroup() != null
                && getReadGroupMaxQual(hdr.getReadGroup(rec.getReadGroup())) != null);
    }

    public static boolean hasFlowTags(final GATKRead rec) {
        return rec.hasAttribute(FlowBasedRead.FLOW_MATRIX_TAG_NAME)
                || rec.hasAttribute(FlowBasedRead.FLOW_MATRiX_OLD_TAG_KR)
                || rec.hasAttribute(FlowBasedRead.FLOW_MATRiX_OLD_TAG_TI);
    }

    public static boolean hasFlowTags(final SAMRecord rec) {
        return rec.hasAttribute(FlowBasedRead.FLOW_MATRIX_TAG_NAME)
                || (rec.getHeader() != null && rec.getAttribute("RG") != null
                    && getReadGroupMaxQual(rec.getHeader().getReadGroup((String)rec.getAttribute("RG"))) != null)
                || rec.hasAttribute(FlowBasedRead.FLOW_MATRiX_OLD_TAG_KR)
                || rec.hasAttribute(FlowBasedRead.FLOW_MATRiX_OLD_TAG_TI);

    }

    /**
     *
     * This is the function to run if you want to ask if the data are flow-based
     *
     * @param hdr - file header
     * @param read - the read
     * @return true if the read is flow-based
     */
    public static boolean isFlowPlatform(final SAMFileHeader hdr, final GATKRead read) {
        if (!hasFlowTags(read, hdr)){
            return false;
        }
        return getReadGroupInfo(hdr, read).isFlowPlatform;
    }
    public static synchronized ReadGroupInfo getReadGroupInfo(final SAMFileHeader hdr, final GATKRead read) {

        if ( !hasFlowTags(read, hdr) ) {
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

    public static int getStrandedUnclippedStartForFlow(final GATKRead read, final SAMFileHeader header, final MarkDuplicatesSparkArgumentCollection mdArgs) {
        return read.isReverseStrand() ? getMarkDupReadEnd(read, false, header, mdArgs) : getMarkDupReadStart(read, false, header, mdArgs);
    }

    public static int getStrandedUnclippedEndForFlow(GATKRead read, SAMFileHeader header, MarkDuplicatesSparkArgumentCollection mdArgs) {
        return !read.isReverseStrand() ? getMarkDupReadEnd(read, true, header, mdArgs) : getMarkDupReadStart(read, true, header, mdArgs);
    }

    /**
     * Get a starting location for the read, mainly for the purpose of comparing it with another one to find duplicates.
     * To reduce uniqueness, there might be different definitions of what a read start is. To begin with, it can be the
     * (hard) clipped or unclipped start. Furthermore, the MarkDuplicates tool defines additional options for determining
     * the start location, such as skipping the first HMER, allowing for uncertainty or the returning of the unclipped
     * location based on a mapping quality value threshold.
     *
     * Note that this function operates in the REFERENCE direction - meaning that for reverse reads it should be called to
     * get the read's end.
     *
     * @param gatkRead - read to get the MarkDuplicates' start location
     * @param endSemantics - location is sought under end-of-fragment semantics
     * @param header - reads file SAMHeader
     * @param mdArgs - MarkDuplicates argument collection
     * @return - read start location, for MarkDuplicates
     */
    public static int getMarkDupReadStart(final GATKRead gatkRead, final boolean endSemantics, final SAMFileHeader header, final MarkDuplicatesSparkArgumentCollection mdArgs) {

        if ( !endSemantics && mdArgs.FLOW_SKIP_START_HOMOPOLYMERS != 0 ) {
            final byte[]      bases = gatkRead.getBasesNoCopy();
            final byte[]      flowOrder = getReadFlowOrder(header, gatkRead);

            byte        hmerBase = bases[0];
            int         flowOrderOfs = 0;
            int         hmersLeft = mdArgs.FLOW_SKIP_START_HOMOPOLYMERS;      // number of hmer left to trim

            // advance flow order to base
            if ( flowOrder != null )
                while ( flowOrder[flowOrderOfs] != hmerBase ) {
                    if (++flowOrderOfs >= flowOrder.length)
                        flowOrderOfs = 0;
                    hmersLeft--;
                }

            int hmerSize = 1;
            for ( ; hmerSize < bases.length ; hmerSize++ )
                if (bases[hmerSize] != hmerBase) {
                    if ( --hmersLeft <= 0 )
                        break;
                    else {
                        hmerBase = bases[hmerSize];
                        if ( flowOrder != null ) {
                            if ( ++flowOrderOfs >= flowOrder.length )
                                flowOrderOfs = 0;
                            while ( flowOrder[flowOrderOfs] != hmerBase ) {
                                hmersLeft--;
                                if ( ++flowOrderOfs >= flowOrder.length )
                                    flowOrderOfs = 0;
                            }
                            if ( hmersLeft <= 0 )
                                break;
                        }
                    }
                }
            final int     start = gatkRead.getUnclippedStart() + hmerSize;
            return mdArgs.FLOW_USE_CLIPPED_LOCATIONS ? Math.max(start, gatkRead.getStart()) : start;
        }
        else if ( readEndMarkedUnclipped(gatkRead, mdArgs.FLOW_Q_IS_KNOWN_END) ) {
            return gatkRead.getUnclippedStart();
        } else if ( endSemantics && readEndMarkedUncertain(gatkRead) ) {
            return FLOW_BASED_INSIGNIFICANT_END;
        } else if ( mdArgs.FLOW_USE_CLIPPED_LOCATIONS ) {
            return gatkRead.getStart();
        } else {
            return gatkRead.getUnclippedStart();
        }
    }

    // this method complements getMarkDupReadStart with respect to the read's end location for MarkDuplicates
    public static int getMarkDupReadEnd(final GATKRead gatkRead, boolean endSemantics, SAMFileHeader header, MarkDuplicatesSparkArgumentCollection mdArgs) {

        if ( !endSemantics && mdArgs.FLOW_SKIP_START_HOMOPOLYMERS != 0 ) {
            final byte[]      bases = gatkRead.getBasesNoCopy();
            final byte[]      flowOrder = getReadFlowOrder(header, gatkRead);

            byte        hmerBase = bases[bases.length - 1];
            int         flowOrderOfs = 0;
            int         hmersLeft = mdArgs.FLOW_SKIP_START_HOMOPOLYMERS;      // number of hmer left to trim

            // advance flow order to base
            if ( flowOrder != null )
                while ( flowOrder[flowOrderOfs] != hmerBase ) {
                    if (++flowOrderOfs >= flowOrder.length)
                        flowOrderOfs = 0;
                    hmersLeft--;
                }

            int         hmerSize = 1;
            for ( ; hmerSize < bases.length ; hmerSize++ )
                if (bases[bases.length - 1 - hmerSize] != hmerBase) {
                    if ( --hmersLeft <= 0 )
                        break;
                    else {
                        hmerBase = bases[bases.length - 1 - hmerSize];
                        if ( flowOrder != null ) {
                            if (++flowOrderOfs >= flowOrder.length)
                                flowOrderOfs = 0;
                            while (flowOrder[flowOrderOfs] != hmerBase) {
                                hmersLeft--;
                                if (++flowOrderOfs >= flowOrder.length)
                                    flowOrderOfs = 0;
                            }
                            if (hmersLeft <= 0)
                                break;
                        }
                    }
                }
            final int     end = gatkRead.getUnclippedEnd() - hmerSize;
            return mdArgs.FLOW_USE_CLIPPED_LOCATIONS ? Math.min(end, gatkRead.getEnd()) : end;
        }
        else if ( readEndMarkedUnclipped(gatkRead, mdArgs.FLOW_Q_IS_KNOWN_END) ) {
            return gatkRead.getUnclippedEnd();
        } else if ( endSemantics && readEndMarkedUncertain(gatkRead) ) {
            return FLOW_BASED_INSIGNIFICANT_END;
        } else if ( mdArgs.FLOW_USE_CLIPPED_LOCATIONS ) {
            return gatkRead.getEnd();
        } else {
            return gatkRead.getUnclippedEnd();
        }
    }

    /**
     * Retrieve flow matrix modifications matrix from its string argument format. This matrix contains
     * logic for modifying the flow matrix as it is read in. If the value of [n] is not zero,
     * then the hmer probability for hmer length n will be copied to the [n] position
     *
     * For the implementation logic, see FlowBasedRead.fillFlowMatrix
     */
    static public int[] getFlowMatrixModsInstructions(final String flowMatrixMods, final int maxHmer) {

        if ( flowMatrixMods != null ) {
            final int[] flowMatrixModsInstructions = new int[maxHmer + 1];

            final String[]    toks = flowMatrixMods.split(",");
            for ( int i = 0 ; i < toks.length - 1 ; i += 2 ) {
                final int hmer = Utils.validIndex(Integer.parseInt(toks[i]), flowMatrixModsInstructions.length);
                flowMatrixModsInstructions[hmer] = Integer.parseInt(toks[i + 1]);
            }
            return flowMatrixModsInstructions;
        } else {
            return null;
        }
    }

    public static CycleSkipStatus getCycleSkipStatus(FlowBasedRead read, ReferenceContext referenceContext) {

        CycleSkipStatus      status = CycleSkipStatus.NS;

        // walk cigar elements of read
        int     readOfs = 0;
        int     refOfs = 0;
        byte[]  bases = read.getBasesNoCopy();
        for (CigarElement cigarElement : read.getCigarElements() ) {
            if ( cigarElement.getOperator() == CigarOperator.M ) {

                // extract bases from element and reference, continue only if different
                byte[]      elemBases = ArrayUtils.subarray(bases, readOfs, readOfs + cigarElement.getLength());
                byte[]      refBases = ArrayUtils.subarray(referenceContext.getBases(), refOfs, refOfs + cigarElement.getLength());
                if ( !Arrays.equals(elemBases, refBases) ) {

                    // create keys
                    int[]  altKey = FlowBasedKeyCodec.baseArrayToKey(elemBases, read.getFlowOrder());
                    int[]  refKey = FlowBasedKeyCodec.baseArrayToKey(refBases, read.getFlowOrder());

                    // assign initial css
                    CycleSkipStatus      cssValue = (refKey.length != altKey.length) ? CycleSkipStatus.CS : CycleSkipStatus.NS;

                    // if same length (NS) then see if it is possible-cycle-skip
                    if ( cssValue == CycleSkipStatus.NS ) {
                        for ( int n = 0 ; n < refKey.length ; n++ ) {
                            if ( (refKey[n] == 0) ^ (altKey[n] == 0) ) {
                                cssValue = CycleSkipStatus.PCS;
                                break;
                            }
                        }
                    }

                    // integrate into status
                    if ( cssValue.getPriority() > status.getPriority() ) {
                        status = cssValue;
                    }
                }

            }
            readOfs += (cigarElement.getOperator().consumesReadBases() ? cigarElement.getLength() : 0);
            refOfs += (cigarElement.getOperator().consumesReferenceBases() ? cigarElement.getLength() : 0);
            if ( status == CycleSkipStatus.CS )
                break;
        }

        return status;
    }

    public static int calcFlowOrderLength(String flowOrder) {

        final int i = flowOrder.indexOf(flowOrder.charAt(0), 1);

        return (i < 0) ? flowOrder.length() : i;
    }
}
