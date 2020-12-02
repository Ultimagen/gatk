package org.ultimagenomics.haplotype_calling;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;

import java.util.*;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

public class LHWRefView {

    final private int             hmerSizeThreshold;

    final private byte[]          fullRef;
    final private Locatable       refLoc;
    final private Logger          logger;
    final private boolean         debug;
    final private SmithWatermanAligner aligner;

    public LHWRefView(final int hmerSizeThreshold, final byte[] fullRef, final Locatable refLoc, final Logger logger, final boolean debug, final SmithWatermanAligner aligner) {

        this.hmerSizeThreshold = hmerSizeThreshold;
        this.fullRef = fullRef;
        this.refLoc = refLoc;
        this.logger = logger;
        this.debug = debug;
        this.aligner = aligner;
        if ( debug ) {
            logger.info("LHWRefView: >" + hmerSizeThreshold + "hmer, refLoc: " + refLoc + " fullRef:");
            logger.info(basesAsString(fullRef));
        }
    }

    public static boolean needsCollapsing(byte[] bases, int hmerSizeThreshold, final Logger logger, final boolean debug) {

        byte    lastBase = 0;
        int     baseSameCount = 0;

        if ( debug ) {
            logger.debug("checking for >" + hmerSizeThreshold + "hmer in:");
            logger.debug(basesAsString(bases));
        }

        // check if has at least one sequence of stable bases larger than threshold
        int         baseIndex = 0;
        for  ( byte base : bases ) {
            if ( base == lastBase ) {
                if ( ++baseSameCount >= hmerSizeThreshold ) {
                    if ( debug )
                        logger.debug("will collapse. found a stable sequence of at least " + (baseSameCount+1) + " of " + Character.toString((char)lastBase));
                    return true;
                }
            } else {
                lastBase = base;
                baseSameCount = 0;
            }
            baseIndex++;
        }

        if ( debug )
            logger.debug("will not collapse");
        return false;
    }

    public List<Haplotype> uncollapseHaplotypesByRef(final Collection<Haplotype> haplotypes, boolean log, boolean limit, byte[] refBases) {

        final List<Haplotype>       result = new LinkedList<>();
        final Map<Locatable, byte[]> refMap = new LinkedHashMap<>();
        int                         alignmentStartHapwrtRef = 0;

        // locate reference haplotype, if needed
        if ( refBases == null )
            for ( Haplotype h : haplotypes )
                if ( h.isReference() ) {
                    refBases = h.getBases();
                    alignmentStartHapwrtRef = h.getAlignmentStartHapwrtRef();
                    break;
                }

        // uncollapse haplotypes
        for ( Haplotype h : haplotypes )
        {
            // by default
            Haplotype   alignedHaplotype = h;

            if ( !h.isReference() ) {

                // find ref for this location
                byte[] ref = refMap.get(h.getGenomeLocation());
                if (ref == null) {
                    ref = uncollapsedPartialRef(h.getGenomeLocation());
                    refMap.put(h.getGenomeLocation(), ref);
                }
                AtomicInteger       offset = new AtomicInteger();
                AtomicBoolean       didCollapse = new AtomicBoolean(false);
                AtomicBoolean       didCollapseRev = new AtomicBoolean(false);
                byte[] alignedBases = uncollapseByRef(h.getBases(), ref, offset, false, didCollapse);
                byte[] alignedBases1 = uncollapseByRef(h.getBases(), ref, offset, true, didCollapseRev);
                if ( alignedBases1.length > alignedBases.length ) {
                    alignedBases = alignedBases1;
                    didCollapse.set(didCollapseRev.get());
                }

                byte[] orgAlignedBases = alignedBases;
                if ( limit )
                    alignedBases = collapseBases(orgAlignedBases);
                alignedHaplotype = new Haplotype(alignedBases, h.isReference());
                alignedHaplotype.setScore(h.getScore());
                alignedHaplotype.setGenomeLocation(h.getGenomeLocation());
                alignedHaplotype.setEventMap(h.getEventMap());
                alignedHaplotype.setAlignmentStartHapwrtRef(offset.get());
                alignedHaplotype.contigs = h.contigs;
                alignedHaplotype.setCollapsed(didCollapse.get());
                alignedHaplotype.setDiffMatter(result.size());
            }

            result.add(alignedHaplotype);
        }

        // if we had a reference, generate cigar against it
        if ( refBases != null ) {
            for ( Haplotype h : result ) {
                if ( !h.isReference() ) {
                    SmithWatermanAlignment  alignment = aligner.align(refBases, h.getBases(), SmithWatermanAligner.ORIGINAL_DEFAULT, SWOverhangStrategy.INDEL);
                    h.setCigar(alignment.getCigar());
                    h.setAlignmentStartHapwrtRef(alignment.getAlignmentOffset() + alignmentStartHapwrtRef);
                }
            }
        }

        if ( log ) {
            int         i = 0;
            for (Haplotype h : haplotypes)
                logHaplotype(h, "COL_" + (limit ? "P1_" : "P2_"), i++);
            i = 0;
            for (Haplotype h : result)
                logHaplotype(h, "UNCOL_" + (limit ? "P1_" : "P2_"), i++);
        }


        return result;
    }

    private byte[] collapseBases(byte[] fullBases) {

        // collapsed sequence would not be longer than full sequence
        byte[]          collapsedBases = new byte[fullBases.length];

        // loop while trimming
        byte    lastBase = 0;
        int     baseSameCount = 0;
        int     srcOfs = 0;
        int     dstOfs = 0;
        boolean firstHomopolymer = true;
        for  ( byte base : fullBases ) {
            if ( base == lastBase ) {
                baseSameCount++;
                if ( !firstHomopolymer && (baseSameCount >= hmerSizeThreshold) ) {
                    // collapsing, do not store
                } else {
                    // stable but under threshold, store
                    collapsedBases[dstOfs++] = base;
                }
            } else {
                // unstable, simply store
                if ( lastBase != 0 )
                    firstHomopolymer = false;
                lastBase = base;
                baseSameCount = 0;

                collapsedBases[dstOfs++] = base;
            }
            srcOfs++;
        }

        // adjust size of collapsedBases
        // not very efficient as it allocates copies the data.
        // do we really need the array to be the right size?
        collapsedBases = Arrays.copyOf(collapsedBases, dstOfs);

        return collapsedBases;
    }

    private byte[] uncollapsedPartialRef(final Locatable ucLoc) {

        int         ucOfs = ucLoc.getStart() - refLoc.getStart();
        int         size = ucLoc.getLengthOnReference();
        byte[]      bases =  Arrays.copyOfRange(fullRef, ucOfs, ucOfs + size);

        if ( debug ) {
            logger.info("uncollapsedPartialRef: cOfs: " + ucOfs + ", size: " + size + ", bases:");
            logger.info(basesAsString(bases));
        }

        return bases;
    }

    private byte[] uncollapseByRef(byte[] bases, byte[] ref, AtomicInteger offsetResult, boolean rev, AtomicBoolean didCollapse) {

        if ( debug ) {
            logger.info("bases, ref, finalResult:");
            logger.info(basesAsString(bases));
            logger.info(basesAsString(ref));
        }

        if ( rev ) {
            byte[]      basesRev = Arrays.copyOf(bases, bases.length);
            byte[]      refRev = Arrays.copyOf(ref, ref.length);
            SequenceUtil.reverseComplement(basesRev);
            SequenceUtil.reverseComplement(refRev);
            bases = basesRev;
            ref = refRev;
        }

        // use aligner to get CIGAR
        SmithWatermanAlignment alignment = aligner.align(ref, bases, SmithWatermanAligner.ORIGINAL_DEFAULT, SWOverhangStrategy.INDEL);
        if ( debug )
            logger.info("alignment.offset: " + alignment.getAlignmentOffset() + ", cigar: " + alignment.getCigar());

        // collect max length by walking the cigar and adding up delete operators (some of which may actually be replaced)
        int     resultLength = bases.length;
        for ( CigarElement c : alignment.getCigar() ) {
            if (c.getOperator() == CigarOperator.D)
                resultLength += c.getLength();
        }
        final byte[] result = new byte[resultLength];

        // prepare offsets
        int         basesOfs = alignment.getAlignmentOffset();
        int         refOfs = 0;
        int         resultOfs = 0;

        // start walking on cigars and make adjustments
        for ( CigarElement c : alignment.getCigar() ) {
            if (c.getOperator() != CigarOperator.D) {
                // not D - simple case
                if (c.getOperator().consumesReadBases()) {
                    System.arraycopy(bases, basesOfs, result, resultOfs, c.getLength());
                    basesOfs += c.getLength();
                    resultOfs += c.getLength();
                }
            } else {

                // check if the incoming bases contain a homopolymer
                int pad = 0;
                int chkPad = 1;
                byte[] fwdSlice = Arrays.copyOfRange(bases, basesOfs, Math.min(basesOfs + hmerSizeThreshold + pad, bases.length));
                byte[] bckSlice = Arrays.copyOfRange(bases, Math.max(0, basesOfs - hmerSizeThreshold - pad), basesOfs);
                if ( needsCollapsing(fwdSlice, hmerSizeThreshold - chkPad, logger, debug) ||
                        needsCollapsing(bckSlice, hmerSizeThreshold - chkPad, logger, debug) ) {


                    // check for a delete at the end of an hmer or at the beginning
                    if (onHomoPolymer(ref, refOfs - hmerSizeThreshold, ref[refOfs], hmerSizeThreshold, 1)) {
                        // fill with base until end of homopolymer on the ref
                        for (int size = 0; (size < c.getLength()) /*&& (ref[refOfs + size] == base)*/; size++)
                            result[resultOfs++] = ref[refOfs + size];
                        didCollapse.set(true);

                    } else if (onHomoPolymer(ref, refOfs + c.getLength(), ref[refOfs + c.getLength() - 1], hmerSizeThreshold, -1)) {
                        // fill with base until start of homopolymer on the ref
                        for (int size = 0; (size < c.getLength()) /*&& (ref[refOfs + c.getLength() - 1 - size] == base)*/; size++)
                            result[resultOfs++] = ref[refOfs + c.getLength() - 1 - size];
                        didCollapse.set(true);
                    }
                }
            }
            if (c.getOperator().consumesReferenceBases())
                refOfs += c.getLength();
        }

        // return adjusted result
        final byte[] finalResult = (result.length == resultOfs) ? result : Arrays.copyOf(result, resultOfs);

        if ( rev )
            SequenceUtil.reverseComplement(finalResult);

        if ( debug ) {
            logger.info(basesAsString(finalResult));
        }

        // return offset?
        if ( offsetResult != null )
            offsetResult.set(alignment.getAlignmentOffset());

        return finalResult;
    }

    private boolean onHomoPolymer(byte[] bases, int ofs, byte base, int length, int direction)
    {
        for ( int tick = 0 ; tick < hmerSizeThreshold ; tick++ ) {
            if ( sameBase(bases, ofs + tick, base, length) )
                return true;
        }
        return false;
    }

    private boolean sameBase(byte[] bases, int ofs, byte base, int length)
    {
        try {
            // has enough bases?
            if (ofs + length > bases.length)
                return false;
            while (length-- != 0) {
                if (bases[ofs++] != base)
                    return false;
            }
            return true;
        } catch (ArrayIndexOutOfBoundsException e) {
            return false;
        }
    }

    private static int logHaplotype_i = 0;
    private void logHaplotype(Haplotype h, String label, int i) {
        if ( debug ) {
            String      name = label + "_" + (logHaplotype_i++) + "_" + i;
            String      contig = h.getGenomeLocation().getContig();
            int         start = h.getGenomeLocation().getStart();
            String      cigar = h.getCigar() != null ? h.getCigar().toString() : "?";
            String      bases = h.getDisplayString();
            if ( h.isReference() )
                name += "_REF";

            byte[]      q = new byte[bases.length()];
            Arrays.fill(q, (byte)40);
            String      quals = new String(q);

            String      line = String.format("%s\t0\t%s\t%d\t60\t%s\t*\t0\t0\t%s\t%s\tRG:Z:%s",
                                                    name, contig, start, cigar, bases, quals, label);

            logger.info(String.format("==SAM==\t%s", line));
        }
    }

    private static String basesAsString(byte[] bases) {
        StringBuilder   sb = new StringBuilder();
        sb.append(String.format("len: %d, bases: ", bases.length));
        sb.append(new String(bases));
        return sb.toString();
    }

}
