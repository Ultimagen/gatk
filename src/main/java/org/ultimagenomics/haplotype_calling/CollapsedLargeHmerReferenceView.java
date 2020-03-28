package org.ultimagenomics.haplotype_calling;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.beanutils.BeanUtils;
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
            logger.info("will not collapse");
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
            //logger.info("collapsedToFullLocationMap: " + Arrays.toString(collapsedToFullLocationMap));
            //logger.info("fullToCollapsedLocationMap: " + Arrays.toString(fullToCollapsedLocationMap));
        }

        collapsedRefLoc = getCollapsedLoc(refLoc);

        // debug: save an alignement beteeen the two references. used for learning and observation. can be removed
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

    public List<VariantContext> uncollapseByRef(List<VariantContext> calls) {

        final List<VariantContext>  result = new LinkedList<>();

        for ( VariantContext other : calls ) {

            // adjust locations
            long        start = toUncollapsedLocus(other.getStart());
            long        end = toUncollapsedLocus(other.getEnd());

            // retrieve true ref content
            int         rangeStart = (int)(start - refLoc.getStart());
            int         rangeSize = other.getEnd() - other.getStart() + 1;
            byte[]      ref = Arrays.copyOfRange(fullRef, rangeStart + 1, rangeStart + rangeSize + 1);

            // modify alleles, getnotype
            List<Allele>        alleles = modifiedAlleles(other.getAlleles(), ref);
            GenotypesContext    genotypes = modifiedGenotypes(other.getGenotypes(), alleles, ref);

            VariantContext      vc = new MyVariantContext(other, start, end, alleles, genotypes);

            result.add(vc);
        }

        for ( VariantContext vc : calls )
            logVariantContext(vc, "uncollapseByRef: >>");
        for ( VariantContext vc : result )
            logVariantContext(vc, "uncollapseByRef: <<");

        return result;
    }

    private List<Allele> modifiedAlleles(List<Allele> otherAlleles, byte[] ref) {

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

    private GenotypesContext modifiedGenotypes(GenotypesContext genotypes, List<Allele> otherAlleles, byte[] ref) {

        GenotypesContext        modGC = GenotypesContext.create();

        for ( Genotype genotype : genotypes ) {

            // modify alleles
            List<Allele>    alleles = modifiedAlleles(genotype.getAlleles(), ref);
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

    public List<Haplotype> uncollapseByRef(final Collection<Haplotype> haplotypes) {

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

            int         loc = h.getGenomeLocation().getStart() - h.getAlignmentStartHapwrtRef();
            int         unloc = toUncollapsedLocus(loc);
            int         unstart = alignedHaplotype.getGenomeLocation().getStart() - unloc;
            alignedHaplotype.setAlignmentStartHapwrtRef(unstart);

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

        for ( Haplotype h : haplotypes )
            logHaplotype(h, "uncollapseByRef: >>");
        for ( Haplotype h : result )
            logHaplotype(h, "uncollapseByRef: <<");


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

    private void logHaplotype(Haplotype h, String label) {
        if ( debug ) {
            logger.info(String.format("%s: %s %d %s %s", label, h.getGenomeLocation(), h.getAlignmentStartHapwrtRef(), h.getCigar(), h));
        }
    }

    private void logVariantContext(VariantContext vc, String label) {
        if ( debug ) {
            logger.info(String.format("%s: %s", label, vc));
        }
    }
}
