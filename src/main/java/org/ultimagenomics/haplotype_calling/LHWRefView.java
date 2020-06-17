package org.ultimagenomics.haplotype_calling;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
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
import java.util.concurrent.atomic.AtomicInteger;

public class LHWRefView {

    private int             hmerSizeThreshold;

    private byte[]          fullRef;
    private Locatable       refLoc;
    private Haplotype       refHaplotype;
    private AssemblyRegion  region;
    private Logger          logger;
    private boolean         debug;
    private byte[]          fullQuals;

    private byte[]          collapsedRef;
    private int[]           fullToCollapsedLocationMap;
    private int[]           collapsedToFullLocationMap;
    private Locatable       collapsedRefLoc;
    private SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);
    SmithWatermanAlignment refAlignement;
    private byte[]          collapsedQuals;

    public LHWRefView(final int hmerSizeThreshold, final byte[] fullRef, final Locatable refLoc, Haplotype refHaplotype, final AssemblyRegion region, final Logger logger, final boolean debug, byte[] fullQuals) {

        this.hmerSizeThreshold = hmerSizeThreshold;
        this.fullRef = fullRef;
        this.refLoc = refLoc;
        this.refHaplotype = refHaplotype;
        this.region = region;
        this.logger = logger;
        this.debug = debug;
        this.fullQuals = fullQuals;

        if ( debug ) {
            logger.info("LHWRefView: >" + hmerSizeThreshold + "hmer, refLoc: " + refLoc + " fullRef:");
            logger.info(printBases(fullRef));
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
        int         baseIndex = 0;
        for  ( byte base : bases ) {
            if ( base == lastBase ) {
                if ( ++baseSameCount >= hmerSizeThreshold ) {
                    if ( debug )
                        logger.info("will collapse. found a stable sequence of at least " + (baseSameCount+1) + " of " + Character.toString((char)lastBase));
                    return true;
                }
            } else {
                lastBase = base;
                baseSameCount = 0;
            }
            baseIndex++;
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

        List<GATKRead>      reads = new LinkedList<>();
        for ( GATKRead read : region.getReads() )
        {
            reads.add(read);
            if (read.getUGOrgLoc() != null )
                continue;

            Locatable       loc = new SimpleInterval(read);
            Locatable       cLoc = getCollapsedLoc(read);
            byte[]          bases = read.getBasesNoCopy();
            boolean         readNeedsCollapsing = needsCollapsing(bases, hmerSizeThreshold, logger, debug);

            if ( debug )
                logger.info(String.format("Read %s %s %c %s -> %s", read.getName(), read.getCigar(),
                                                    readNeedsCollapsing ? 'C' : 'N',
                                                    loc, cLoc));
            if ( readNeedsCollapsing ) {
                byte[]          quals = read.getBaseQualitiesNoCopy();
                LHWRefView      assist = new LHWRefView(hmerSizeThreshold, read.getBases(), null, null, null, logger, debug, quals);
                byte[]          cBases = assist.collapsedRef;
                byte[]          cQuals = assist.collapsedQuals;

                if ( debug ) {
                    logger.info(printBases(bases));
                    logger.info(printBases(quals));
                    logger.info(printBases(cBases));
                    logger.info(printBases(cQuals));
                }

                // DK!!! Big Hack!
                // instead of adjusting the Cigar, the number of bases is kept, leaving a wrong cigar, but
                // at least with the right length.
                // Also - Qualities seems to be always 40 - why?
                byte[]          cBases1 = extendByteArray(cBases, bases.length, (byte)'N');
                byte[]          cQuals1 = extendByteArray(cQuals, quals.length, (byte)40);
                read.setUGOrgBases(bases);
                read.setUGOrgQuals(quals);
                read.setBases(cBases1);
                read.setBaseQualities(cQuals1);
            }
            read.setUGOrgPosition(loc);
            read.setPosition(cLoc);
        }

        return reads;
    }

    private byte[] extendByteArray(byte[] src, int len, byte fillValue) {

        byte[]          dst = new byte[len];
        System.arraycopy(src, 0, dst, 0, src.length);
        Arrays.fill(dst, src.length, len, fillValue);

        return dst;
    }

    private void collapse() {

        // collapsed sequence would not be longer than full sequence
        collapsedRef = new byte[fullRef.length];
        fullToCollapsedLocationMap = new int[fullRef.length];
        collapsedToFullLocationMap = new int[fullRef.length];
        if ( fullQuals != null )
            collapsedQuals = new byte[fullRef.length];
        else
            collapsedQuals = null;

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
                    if ( fullQuals != null )
                        collapsedQuals[dstOfs] = fullQuals[srcOfs];
                    collapsedRef[dstOfs++] = base;
                }
            } else {
                // unstable, simply store
                lastBase = base;
                baseSameCount = 0;
                fullToCollapsedLocationMap[srcOfs] = dstOfs;
                collapsedToFullLocationMap[dstOfs] = srcOfs;
                if ( fullQuals != null )
                    collapsedQuals[dstOfs] = fullQuals[srcOfs];

                // just coming out of collapsing?
                if ( dstOfs > 0 && ((collapsedToFullLocationMap[dstOfs - 1] + 1) != collapsedToFullLocationMap[dstOfs]) ) {
                    // attribute the second half of the collapsed area to end of the original area
                    for ( int n = 1 ; n <= (hmerSizeThreshold / 2) ; n++ )
                        collapsedToFullLocationMap[dstOfs - n] = collapsedToFullLocationMap[dstOfs] - n;
                }

                collapsedRef[dstOfs++] = base;
            }
            srcOfs++;
        }

        // adjust size of collapsedRef
        // not very efficient as it allocates copies the data.
        // do we really need the array to be the right size?
        collapsedRef = Arrays.copyOf(collapsedRef, dstOfs);
        collapsedToFullLocationMap = Arrays.copyOf(collapsedToFullLocationMap, dstOfs);
        if ( collapsedQuals != null )
            collapsedQuals = Arrays.copyOf(collapsedQuals, dstOfs);

        if ( debug ) {
            logger.info("after collapsing: ");
            logger.info(printBases(collapsedRef));
            //logger.info("collapsedToFullLocationMap: " + Arrays.toString(collapsedToFullLocationMap));
            //logger.info("fullToCollapsedLocationMap: " + Arrays.toString(fullToCollapsedLocationMap));
        }

        if ( refLoc != null )
            collapsedRefLoc = getCollapsedLoc(refLoc);
    }

    private byte[] collapseBases(byte[] fullBases) {

        // collapsed sequence would not be longer than full sequence
        byte[]          collapsedBases = new byte[fullBases.length];

        // loop while trimming
        byte    lastBase = 0;
        int     baseSameCount = 0;
        int     srcOfs = 0;
        int     dstOfs = 0;
        for  ( byte base : fullBases ) {
            if ( base == lastBase ) {
                if ( ++baseSameCount >= hmerSizeThreshold ) {
                    // collapsing, do not store
                } else {
                    // stable but under threshold, store
                    collapsedBases[dstOfs++] = base;
                }
            } else {
                // unstable, simply store
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

    public byte[] uncollapseByRef(byte[] bases, byte[] ref, AtomicInteger offsetResult) {

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

        // return offset?
        if ( offsetResult != null )
            offsetResult.set(alignment.getAlignmentOffset());

        return finalResult;
    }

    public List<VariantContext> uncollapseCallsByRef(List<VariantContext> calls) {

        final List<VariantContext>  result = new LinkedList<>();

        for ( VariantContext other : calls ) {

            // adjust locations
            long        start = toUncollapsedLocus(other.getStart());
            long        end = toUncollapsedLocus(other.getEnd());

            // retrieve true ref content
            int rangeStart = (int) (start - refLoc.getStart());
            int rangeSize = other.getEnd() - other.getStart() + 1;
            byte[] ref = Arrays.copyOfRange(fullRef, rangeStart + 1, rangeStart + rangeSize + 1);

            // modify alleles, getnotype
            List<Allele> alleles = replaceRefInAlleles(other.getAlleles(), ref);
            GenotypesContext genotypes = replaceRefInGenotypes(other.getGenotypes(), ref);

            VariantContext vc = new MyVariantContext(other, start, end, alleles, genotypes);

            result.add(vc);
        }

        for ( VariantContext vc : calls )
            logVariantContext(vc, "uncollapseByRef: >>");
        for ( VariantContext vc : result )
            logVariantContext(vc, "uncollapseByRef: <<");

        return result;
    }

    public List<Haplotype> uncollapseHaplotypesByRef(final Collection<Haplotype> haplotypes, boolean log, boolean limit) {

        final List<Haplotype>       result = new LinkedList<>();
        final Map<Locatable, byte[]> refMap = new LinkedHashMap<>();
        Haplotype                   refHaplotype = null;

        // locate reference haplotype, if exists
        for ( Haplotype h : haplotypes )
            if ( h.isReference() ) {
                refHaplotype = h;
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
                    ref = getUncollapsedPartialRef(h.getGenomeLocation(), true);
                    refMap.put(h.getGenomeLocation(), ref);
                }
                AtomicInteger       offset = new AtomicInteger();
                byte[] alignedBases = uncollapseByRef(h.getBases(), ref, offset);
                if ( limit )
                    alignedBases = collapseBases(alignedBases);
                alignedHaplotype = new Haplotype(alignedBases, h.isReference());
                alignedHaplotype.setScore(h.getScore());
                alignedHaplotype.setGenomeLocation(getUncollapsedLoc(h.getGenomeLocation()));
                alignedHaplotype.setEventMap(h.getEventMap());
                alignedHaplotype.setAlignmentStartHapwrtRef(offset.get());
            }

            result.add(alignedHaplotype);
        }

        // if we had a reference, generate cigar against it
        if ( refHaplotype != null ) {
            for ( Haplotype h : result ) {
                if ( h != refHaplotype ) {
                    SmithWatermanAlignment  alignment = aligner.align(refHaplotype.getBases(), h.getBases(), SmithWatermanAligner.ORIGINAL_DEFAULT, SWOverhangStrategy.INDEL);
                    h.setCigar(alignment.getCigar());
                    h.setAlignmentStartHapwrtRef(alignment.getAlignmentOffset() + refHaplotype.getAlignmentStartHapwrtRef());
                }
            }
        }

        if ( log ) {
            for (Haplotype h : haplotypes)
                logHaplotype(h, "COL");
            for (Haplotype h : result)
                logHaplotype(h, "UNCOL");
        }


        return result;
    }


    private List<Allele> replaceRefInAlleles(List<Allele> otherAlleles, byte[] ref) {

        // modify alleles
        List<Allele>    alleles = new LinkedList<>();
        Allele          refAllele = null;
        for ( Allele otherAllele : otherAlleles ) {
            if (otherAllele.isReference()) {

                refAllele = Allele.create(ref, true);
                alleles.add(refAllele);

            } else
                alleles.add(otherAllele);
        }

        return alleles;
    }

    private GenotypesContext replaceRefInGenotypes(GenotypesContext genotypes, byte[] ref) {

        GenotypesContext        modGC = GenotypesContext.create();

        for ( Genotype genotype : genotypes ) {

            // modify alleles
            List<Allele>    alleles = replaceRefInAlleles(genotype.getAlleles(), ref);
            Genotype        g = GenotypeBuilder.create(genotype.getSampleName(), alleles, genotype.getExtendedAttributes());

            modGC.add(g);
        }

        return modGC;
    }

    public static class MyVariantContext extends VariantContext {

        public static final long serialVersionUID = 1L;

        public MyVariantContext(VariantContext other, long start, long end, List<Allele> alleles, GenotypesContext genotypes) {

            super(other.getSource(), other.getID(), other.getContig(), start, end, alleles, genotypes, other.getLog10PError(), other.getFiltersMaybeNull(), other.getAttributes(),
                    /*other.fullyDecoded*/ false,
                    EnumSet.noneOf(VariantContext.Validation.class));
        }
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

    private static int logHaplotype_i = 0;
    private void logHaplotype(Haplotype h, String label) {
        if ( debug ) {
            String      name = label + "_" + (logHaplotype_i++);
            String      contig = h.getGenomeLocation().getContig();
            int         start = h.getGenomeLocation().getStart();
            String      cigar = h.getCigar() != null ? h.getCigar().toString() : "?";
            String      bases = h.getDisplayString();

            byte[]      q = new byte[bases.length()];
            Arrays.fill(q, (byte)40);
            String      quals = new String(q);

            String      line = String.format("%s\t0\t%s\t%d\t60\t%s\t*\t0\t0\t%s\t%s\tRG:Z:%s",
                                                    name, contig, start, cigar, bases, quals, label);

            logger.info(String.format("==SAM==\t%s", line));
        }
    }

    private void logVariantContext(VariantContext vc, String label) {
        if ( debug ) {
            logger.info(String.format("%s: %s", label, vc));
        }
    }

    public byte[] getFullRef() {
        return fullRef;
    }

    public Locatable getRefLoc() {
        return refLoc;
    }

}
