package org.ultimagenomics.haplotype_calling;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;

import java.util.*;
import java.util.stream.StreamSupport;

public class CollapsedLargeHmerReferenceView {

    private int             hmerSizeThreshold;
    private byte[]          fullRef;
    private Locatable       refLoc;
    private Haplotype       refHaplotype;
    private AssemblyRegion  region;
    private Logger          logger;
    private boolean         debug;

    private byte[]          collapsedRef;
    private int[]           fullToCollapsedLocationMap;
    private int[]           collapsedToFullLocationMap;
    private Locatable       collapsedRefLoc;
    private SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);
    SmithWatermanAlignment refAlignement;

    public CollapsedLargeHmerReferenceView(final int hmerSizeThreshold, final byte[] fullRef, final Locatable refLoc, Haplotype refHaplotype, final AssemblyRegion region, final Logger logger, final boolean debug) {

        this.hmerSizeThreshold = hmerSizeThreshold;
        this.fullRef = fullRef;
        this.refLoc = refLoc;
        this.refHaplotype = refHaplotype;
        this.region = region;
        this.logger = logger;
        this.debug = debug;

        if ( debug ) {
            logger.info("CollapsedLargeHmerReferenceView: >" + hmerSizeThreshold + "hmer, refLoc: " + refLoc + " fullRef:");
            logger.info(printBases(fullRef));
        }

        // TEMP! test if actual refHaplotype needs collapsing to aid in finding interesting areas for testing
        if ( needsCollapsing(refHaplotype.getBases(), hmerSizeThreshold, logger, debug) ) {
            if ( debug )
                logger.info("refHaplotype needs collapsing!");
        }

        collapse();
    }

    public static boolean needsCollapsing(byte[] bases, int hmerSizeThreshold, final Logger logger, final boolean debug) {

        byte    lastBase = 0;
        int     baseSameCount = 0;

        if ( debug ) {
            logger.info("checking for >" + hmerSizeThreshold + "hmer in:");
            logger.info(printBases(bases));
        }

        // check if has at least one sequence of stable bases larger than threshold
        for  ( byte base : bases ) {
            if ( base == lastBase ) {
                if ( ++baseSameCount > hmerSizeThreshold ) {
                    if ( debug )
                        logger.info("will collapse. found a stable sequence of at least " + baseSameCount + " of " + Character.toString((char)lastBase));
                    return true;
                }
            } else {
                lastBase = base;
                baseSameCount = 0;
            }
        }

        if ( debug )
            logger.info("will not ollapse");
        return false;
    }

    public byte[] getCollapsedFullRef() {
        Utils.nonNull(collapsedRef);
        return collapsedRef;
    }

    public SimpleInterval getCollapsedLoc(final Locatable loc) {
        Utils.nonNull(collapsedRef);
        SimpleInterval  cLoc =  new SimpleInterval(loc.getContig(), toCollapsedLocus(loc.getStart()), toCollapsedLocus(loc.getEnd()));

        if ( debug )
            logger.info("getCollapsedLoc: " + loc + " -> " + cLoc);

        return cLoc;
    }

    public SimpleInterval getUncollapsedLoc(final Locatable loc) {
        Utils.nonNull(collapsedRef);
        SimpleInterval  ucLoc =  new SimpleInterval(loc.getContig(), toUncollapsedLocus(loc.getStart()), toUncollapsedLocus(loc.getEnd()));

        if ( debug )
            logger.info("getUncollapsedLoc: " + loc + " -> " + ucLoc);

        return ucLoc;
    }

    public byte[] getCollapsedPartialRef(final Locatable loc) {
        Utils.nonNull(collapsedRef);

        Locatable   cLoc = getCollapsedLoc(loc);
        int         cOfs = cLoc.getStart() - refLoc.getStart();
        int         size = cLoc.getLengthOnReference();
        byte[]      bases =  Arrays.copyOfRange(collapsedRef, cOfs, cOfs + size);

        if ( debug ) {
            logger.info("getCollapsedPartialRef: cOfs: " + cOfs + ", size: " + size + ", bases:");
            logger.info(printBases(bases));
        }

        return bases;
    }

    public byte[] getUncollapsedPartialRef(final Locatable loc, boolean uncollapseLoc) {
        Utils.nonNull(collapsedRef);

        Locatable   ucLoc = uncollapseLoc ? getUncollapsedLoc(loc) : loc;
        int         ucOfs = ucLoc.getStart() - refLoc.getStart();
        int         size = ucLoc.getLengthOnReference();
        byte[]      bases =  Arrays.copyOfRange(fullRef, ucOfs, ucOfs + size);

        if ( debug ) {
            logger.info("getUncollapsedPartialRef: cOfs: " + ucOfs + ", size: " + size + ", bases:");
            logger.info(printBases(bases));
        }

        return bases;
    }

    public Haplotype getCollapsedRefHaplotype(final Locatable loc) {
        Utils.nonNull(collapsedRef);

        // calc start/end and bases in collapsed space
        Locatable       activeRegion = getCollapsedLoc(loc);
        byte[]          refBases = getCollapsedPartialRef(loc);

        final Haplotype viewHaplotype = new Haplotype(refBases, true);
        viewHaplotype.setAlignmentStartHapwrtRef(activeRegion.getStart() - refLoc.getStart());
        final Cigar c = new Cigar();
        c.add(new CigarElement(viewHaplotype.getBases().length, CigarOperator.M));
        viewHaplotype.setCigar(c);

        if ( debug ) {
            logger.info("getCollapsedRefHaplotype: viewHaplotype: " + viewHaplotype);
            logger.info("getCollapsedRefHaplotype: refHaplotype:  " + refHaplotype);
        }
        return viewHaplotype;
    }

    static final boolean HARD_COLLAPSE_READS = true;

    public List<GATKRead> getCollapsedReads(final AssemblyRegion region) {
        Utils.nonNull(collapsedRef);

        final List<GATKRead>      reads = region.getReads();
        for ( GATKRead read : reads )
            read.setPosition(getCollapsedLoc(read));
        return reads;
    }

    private void collapse() {

        // collapsed sequence would not be longer than full sequence
        collapsedRef = new byte[fullRef.length];
        fullToCollapsedLocationMap = new int[fullRef.length];
        collapsedToFullLocationMap = new int[fullRef.length];

        // loop while trimming
        byte    lastBase = 0;
        int     baseSameCount = 0;
        int     srcOfs = 0;
        int     dstOfs = 0;
        for  ( byte base : fullRef ) {
            if ( base == lastBase ) {
                if ( ++baseSameCount >= hmerSizeThreshold ) {
                    // collapsing, do not store
                    fullToCollapsedLocationMap[srcOfs] = dstOfs - 1;
                } else {
                    // stable but under threshold, store
                    fullToCollapsedLocationMap[srcOfs] = dstOfs;
                    collapsedToFullLocationMap[dstOfs] = srcOfs;
                    collapsedRef[dstOfs++] = base;
                }
            } else {
                // unstable, simply store
                lastBase = base;
                baseSameCount = 0;
                fullToCollapsedLocationMap[srcOfs] = dstOfs;
                collapsedToFullLocationMap[dstOfs] = srcOfs;
                collapsedRef[dstOfs++] = base;
            }
            srcOfs++;
        }

        // adjust size of collapsedRef
        // not very efficient as it allocates copies the data.
        // do we really need the array to be the right size?
        collapsedRef = Arrays.copyOf(collapsedRef, dstOfs);
        collapsedToFullLocationMap = Arrays.copyOf(collapsedToFullLocationMap, dstOfs);

        if ( debug ) {
            logger.info("after collapsing: ");
            logger.info(printBases(collapsedRef));
            logger.info("collapsedToFullLocationMap: " + Arrays.toString(collapsedToFullLocationMap));
            logger.info("fullToCollapsedLocationMap: " + Arrays.toString(fullToCollapsedLocationMap));
        }

        collapsedRefLoc = getCollapsedLoc(refLoc);

        // debug: save an alignement beteeen the two references. used for learnign and observation. can be removed
        refAlignement = aligner.align(fullRef, collapsedRef, SmithWatermanAligner.ORIGINAL_DEFAULT, SWOverhangStrategy.INDEL);
        if ( debug ) {
            logger.info("alignment.offset: " + refAlignement.getAlignmentOffset() + ", cigar: " + refAlignement.getCigar());
        }
        uncollapseByRef(collapsedRef, fullRef);
    }
    
    private int toUncollapsedLocus(int locus) {
        return refLoc.getStart() + collapsedToFullLocationMap[locus - refLoc.getStart()];
    }

    private int toCollapsedLocus(int locus) {
        return refLoc.getStart() + fullToCollapsedLocationMap[locus - refLoc.getStart()];
    }

    private static String printBases(byte[] bases) {
        StringBuilder   sb = new StringBuilder();
        sb.append(String.format("len: %d, bases: ", bases.length));
        sb.append(new String(bases));
        return sb.toString();
    }

    public byte[] uncollapseByRef(byte[] bases, byte[] ref) {

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
                // on ref, check if D has atleast threshold same hmers to the left and D.length to the right
                // if so, this is a candidate for uncollapsing
                byte        base = ref[refOfs];
                if ( sameBase(ref, refOfs, base, hmerSizeThreshold + c.getLength()) ) {
                    Arrays.fill(result, resultOfs, resultOfs + c.getLength(), base);
                    resultOfs += c.getLength();
                }
            }
            if (c.getOperator().consumesReferenceBases())
                refOfs += c.getLength();
        }

        // return adjusted result
        final byte[] finalResult = (result.length == resultOfs) ? result : Arrays.copyOf(result, resultOfs);

        if ( debug ) {
            logger.info("bases, ref, finalResult:");
            logger.info(printBases(bases));
            logger.info(printBases(ref));
            logger.info(printBases(finalResult));
        }

        return finalResult;
    }

    public List<Haplotype> uncollapseByRef(final List<Haplotype> haplotypes) {

        final List<Haplotype>       result = new LinkedList<>();
        final Map<Locatable, byte[]> refMap = new LinkedHashMap<>();
        Haplotype                   refHaplotype = null;

        // uncollapse haplotypes
        for ( Haplotype h : haplotypes )
        {
            // find ref for this location
            byte[]      ref = refMap.get(h.getGenomeLocation());
            if ( ref == null ) {
                ref = getUncollapsedPartialRef(h.getGenomeLocation(), true);
                refMap.put(h.getGenomeLocation(), ref);
            }
            byte[]      alignedBases = uncollapseByRef(h.getBases(), ref);
            Haplotype   alignedHaplotype = new Haplotype(alignedBases, h.isReference());
            alignedHaplotype.setScore(h.getScore());
            alignedHaplotype.setGenomeLocation(getUncollapsedLoc(h.getGenomeLocation()));
            alignedHaplotype.setEventMap(h.getEventMap());

            // save refHaplotype?
            if ( refHaplotype == null && alignedHaplotype.isReference() )
                refHaplotype = alignedHaplotype;

            result.add(alignedHaplotype);
        }

        // if we had a reference, generate cigar against it
        if ( refHaplotype != null ) {
            for ( Haplotype h : result ) {

                SmithWatermanAlignment  alignment = aligner.align(refHaplotype.getBases(), h.getBases(), SmithWatermanAligner.ORIGINAL_DEFAULT, SWOverhangStrategy.INDEL);
                h.setCigar(alignment.getCigar());
            }
        }

        return result;
    }


    private boolean sameBase(byte[] bases, int ofs, byte base, int length)
    {
        // has enough bases?
        if ( ofs + length > bases.length )
            return false;
        while ( length-- != 0 ) {
            if (bases[ofs++] != base)
                return false;
        }
        return true;
    }
}
