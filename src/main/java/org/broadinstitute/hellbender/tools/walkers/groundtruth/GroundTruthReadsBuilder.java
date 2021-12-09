package org.broadinstitute.hellbender.tools.walkers.groundtruth;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.Tuple;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.ToolSuccessfulPrematureExit;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.tools.walkers.featuremapping.FlowFeatureMapper;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.FlowBasedAlignmentEngine;
import org.broadinstitute.hellbender.utils.haplotype.FlowBasedHaplotype;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import picard.cmdline.programgroups.BaseCallingProgramGroup;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Random;
import java.util.zip.GZIPOutputStream;

/**
 * An internal tool to produce a flexible and robust ground truth set for base calling training.
 *
 *
 * <h3> Input </h3>
 * <ul>
 *     <li> Coordinate-sorted and indexed SAM/BAM/CRAM </li>
 *     <li> Maternal and Parental references (fa) </li>
 *     <li> Folder with address translation files from reference to maternal/parental references (filename example: maternal.chr9.csv)</li>
 * </ul>
 *
 * <h3> Output </h3>
 * <ul>
 *     <li> CSV file containing maternal/parental haplotype scores and many more columns </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 * <pre>
 * gatk GroundTruthReadsBuilder \
 *  -R
 *  ../../../ref/Homo_sapiens_assembly38.fasta
 *  -I
 *  150548-UGAv3-4.chr9.cram
 *  --maternal-ref
 *  chr9_HG001_maternal.fa
 *  --paternal-ref
 *  chr9_HG001_paternal.fa
 *  --ancestral-translators-base-path
 *  ./
 *  --output-csv
 *  output-small.csv
 *  --subsampling-ratio
 *  1.0
 *  --max-output-reads
 *  100000000
 *  --intervals
 *  chr9:109991494-109991494
 *  --smith-waterman
 *  FASTEST_AVAILABLE
 *  --likelihood-calculation-engine
 *  FlowBased
 *  -mbq
 *  0
 *  --kmer-size
 *  10
 *  --gt-debug
 *  --output-flow-length
 *  1000
 *  --haplotype-output-padding-size
 *  8
 *  --prepend-sequence
 *  TTTT
 *  --append-sequence
 *  CCCC
 * </pre>
 *
 * {@GATK.walkertype ReadWalker}
 */
@CommandLineProgramProperties(
        summary = "Ground Truth Reads Builder",
        oneLineSummary = "Pproduces a flexible and robust ground truth set for base calling training",
        programGroup = BaseCallingProgramGroup.class
)

@DocumentedFeature
@ExperimentalFeature
public final class GroundTruthReadsBuilder extends ReadWalker {

    private static final Logger logger = LogManager.getLogger(GroundTruthReadsBuilder.class);
    public static final int DEFAULT_FILL_VALUE = -65;
    public static final int NONREF_FILL_VALUE = -80;
    public static final int UNKNOWN_FILL_VALUE = -85;
    public static final int SOFTCLIP_FILL_VALUE = -83;
    private static final int EXTRA_FILL_FROM_HAPLOTYPE = 50;

    @ArgumentCollection
    private final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();

    @Argument(fullName = "maternal-ref", doc="maternal reference file")
    public GATKPath maternalRefPath = null;
    @Argument(fullName = "paternal-ref", doc="paternal reference file")
    public GATKPath paternalRefPath = null;
    @Argument(fullName = "ancestral-translators-base-path", doc="base path for ancestral translation ancestral.contig.csv files")
    public GATKPath ancestralTranslatorsBasePath = null;

    @Argument(fullName = "subsampling-ratio", doc = "subsampling ratio, should be between 0 and 1", optional = true)
    public double   subsamplingRatio = 1.0;
    @Argument(fullName = "max-output-reads", doc = "maximal number of reads to output", optional = true)
    public int      maxOutputReads = 20000000;

    @Argument(fullName = "output-flow-length", doc = "Required length of output flows", optional = true)
    public int      outputFlowLength = 0;
    @Argument(fullName = "prepend-sequence", doc = "Sequence to prepend (barcode)", optional = true)
    public String   prependSequence;
    @Argument(fullName = "append-sequence", doc = "Sequence to append (adapter)", optional = true)
    public String   appendSequence;

    @Argument(fullName = "min-mq", doc = "Minimal mapping quality", optional = true)
    public double   minMappingQuality = 0;
    @Argument(fullName = "max-rq", doc = "Maximal read quality", optional = true)
    public double   maxReadQuality = 0;
    @Argument(fullName = "include-supp-align", doc = "Include supplementary alignments", optional = true)
    public boolean  includeSuppAlign = false;
    @Argument(fullName = "min-haplotype-score", doc = "Minimal score (likelihood) on either haplotype", optional = true)
    public double   minHaplotypeScore = 0;
    @Argument(fullName = "min-haplotype-score-delta", doc = "Minimal score (likelihood) delta between haplotypes", optional = true)
    public double   minHaplotypeScoreDelta = 0;
    @Argument(fullName = "haplotype-output-padding-size", doc = "Number of N to append to best haplotype on output", optional = true)
    public int      haplotypeOutputPaddingSize = 8;
    @Argument(fullName = "discard-non-polyt-softclipped-reads", doc = "Discard reads which are softclipped, unless the softclip is polyT, defaults to true", optional = true)
    public boolean  discardNonPolytSoftclippedReads = false;

    @Argument(fullName = "fill-trimmed-reads-Q", doc = "Reads with tm:Q should be filled from haplotype, otherwise (default) filled with -80", optional = true)
    public boolean fillTrimmedReadsQ;
    @Argument(fullName = "fill-trimmed-reads-Z", doc = "Reads with tm:Z should be filled from haplotype, otherwise (default) filled with -80", optional = true)
    public boolean fillTrimmedReadsZ;
    @Argument(fullName = "fill-trimmed-reads", doc = "Reads with tm:Q or tm:Z should be filled from haplotype, otherwise (default) filled with -80", optional = true)
    public boolean fillTrimmedReads;
    @Argument(fullName = "fill-softclipped-reads", doc = "Softclipped reads should be filled from haplotype, otherwise (default) filled with -83", optional = true)
    public boolean fillSoftclippedReads;

    @Argument(fullName = "output-csv", doc="main output file")
    public GATKPath outputCsvPath = null;

    @Argument(fullName = "gt-debug", doc = "Turn additional internal logging on", optional = true)
    public boolean      debugMode = false;

    @Argument(fullName = "gt-no-output", doc = "do not generate output records", optional = true)
    public boolean      noOutput = false;

    // locals
    final Random                        random = new Random();
    int                                 outputReadsCount = 0;
    ReferenceDataSource                 maternalReference;
    ReferenceDataSource                 paternalReference;
    AncestralContigLocationTranslator   locationTranslator;
    FlowBasedAlignmentEngine            likelihoodCalculationEngine;
    PrintWriter                         outputCsv;
    private final Map<String, ReadGroupInfo> readGroupInfo = new LinkedHashMap<>();
    private int                         locationTranslationErrors;

    // static/const
    static final private String[]       CSV_FIELD_ORDER = {
            "ReadName", "ReadChrom", "ReadStart", "ReadEnd",
            "PaternalHaplotypeScore", "MaternalHaplotypeScore", "RefHaplotypeScore",
            "ReadKey", "BestHaplotypeKey", "ConsensusHaplotypeKey",
            "tm", "mapq", "flags", "ReadCigar",
            "ReadSequence", "PaternalHaplotypeSequence", "MaternalHaplotypeSequence", "BestHaplotypeSequence",
            "ReadUnclippedStart", "ReadUnclippedEnd", "PaternalHaplotypeInterval", "MaternalHaplotypeInterval"
    };


    private static class ScoredHaplotype {
        ReferenceContext        ref;
        ReferenceContext        clippedRef;
        ReferenceContext        unclippedRef;
        int                     softclipFrontFillCount;
        Haplotype               haplotype;
        double                  score;
    }

    static private class ReadGroupInfo {
        public String  flowOrder;
        public String  reversedFlowOrder;
        public int     maxClass;
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        // initialize references
        maternalReference = ReferenceDataSource.of(maternalRefPath.toPath());
        paternalReference = ReferenceDataSource.of(paternalRefPath.toPath());
        locationTranslator = new AncestralContigLocationTranslator(ancestralTranslatorsBasePath);

        // create likelihood engine
        ReadLikelihoodCalculationEngine engine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(hcArgs.likelihoodArgs, false);
        if ( engine instanceof FlowBasedAlignmentEngine) {
            likelihoodCalculationEngine = (FlowBasedAlignmentEngine)engine;
        } else {
            throw new GATKException("must use a flow based likelihood calculation engine");
        }

        // open output, write header
        try {
            if (outputCsvPath.toPath().toString().endsWith(".gz")) {
                outputCsv = new PrintWriter(new GZIPOutputStream(outputCsvPath.getOutputStream()));
            } else {
                outputCsv = new PrintWriter(outputCsvPath.getOutputStream());
            }
        } catch (IOException e) {
            throw new GATKException("failed to open csv output: " + outputCsvPath, e);
        }
        emitCsvHeaders();
    }

    @Override
    public void closeTool() {

        if ( locationTranslationErrors != 0 ) {
            logger.warn("" + locationTranslationErrors + " location translation errors detected");
        }

        outputCsv.close();
        super.closeTool();
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        // filter out due to mapping quality
        if ( minMappingQuality != 0 && read.getMappingQuality() < minMappingQuality ) {
            return;
        }

        // supplemental alignment filter
        if ( read.isSupplementaryAlignment() && !includeSuppAlign ) {
            return;
        }

        // discard because softclipped
        if ( discardNonPolytSoftclippedReads && isEndSoftclipped(read) && !isEndPolyTSoftclipped(read) ) {
            return;
        }

        // subsample
        // NOTE: this is done BEFORE read quality and haplotype scoring
        if ( random.nextDouble() > subsamplingRatio ) {
            return;
        }

        // filter out due to read quality
        FlowBasedRead                   flowRead = null;
        if ( maxReadQuality != 0 ) {
            flowRead = buildFlowRead(read);
            if (getFlowBasedReadQuality(flowRead, flowRead.getMaxHmer()) > maxReadQuality) {
                return;
            }
        }

        // process the read
        try {
            // make sure we have a flow read
            if ( flowRead == null ) {
                flowRead = buildFlowRead(read);
            }

            // prepare
            final ScoredHaplotype                   maternal = new ScoredHaplotype();
            final ScoredHaplotype                   paternal = new ScoredHaplotype();

            // translate location to ascentors
            final Tuple<SimpleInterval, SimpleInterval>  ancestralLocs = locationTranslator.translate(read);
            maternal.ref = new ReferenceContext(maternalReference, ancestralLocs.a);
            paternal.ref = new ReferenceContext(paternalReference, ancestralLocs.b);

            // build haplotypes
            maternal.haplotype = buildReferenceHaplotype(maternal.ref, read);
            paternal.haplotype = buildReferenceHaplotype(paternal.ref, read);
            buildExtendedRef(maternal, maternalReference, ancestralLocs.a, read);
            buildExtendedRef(paternal, paternalReference, ancestralLocs.b, read);

            // generate score for reference
            final double        refScore = scoreReadAgainstReference(read, referenceContext);

            // score read against haplotypes, create flow versions of read nad haplotype
            if ( areSame(maternal.haplotype, referenceContext, read.isReverseStrand()) ) {
                maternal.score = refScore;
            } else {
                maternal.score = scoreReadAgainstHaplotype(read, maternal);
            }
            if ( areSame(paternal.haplotype, referenceContext, read.isReverseStrand()) ) {
                paternal.score = refScore;
            } else if ( arsSame(maternal.haplotype, paternal.haplotype, read.isReverseStrand()) ) {
                paternal.score = scoreReadAgainstHaplotype(read, paternal);
            } else {
                paternal.score = scoreReadAgainstHaplotype(read, paternal);
            }

            // debug printing (in INFO for now, will be changed to DEBUG)
            debugLog(read, referenceContext, maternal, paternal);

            // filter on min score
            // TODO: this is probaby wrong since the scores are negative. To be handled later
            if ( minHaplotypeScore != 0 && Math.min(maternal.score, paternal.score) > minHaplotypeScore ) {
                return;
            }

            // filter on score delta
            if ( minHaplotypeScoreDelta != 0 && Math.abs(maternal.score - paternal.score) > minHaplotypeScoreDelta ) {
                return;
            }

            // if here, emit this read
            outputReadsCount++;
            emit(read, flowRead, refScore, maternal, paternal);

            // limit number of output reads
            if ( maxOutputReads != 0 && outputReadsCount >= maxOutputReads ) {
                // terminate tool
                throw new ToolSuccessfulPrematureExit("stopping tool, output reads max reached: \" + maxOutputReads");
            }

        } catch (LocationTranslationException e) {
            logger.warn("location translation exception: " + e.getMessage());
            locationTranslationErrors++;
        } catch (IOException e) {
            throw new GATKException("failed to process read: " + read.getName(), e);
        }
    }

    private boolean shouldFillFromHaplotype(final GATKRead read) {

        // softclip has priori
        if ( isEndSoftclipped(read) )
            return fillSoftclippedReads;

        // extending timmed as well?
        final String    tm = read.getAttributeAsString("tm");
        if ( tm == null ) {
            return true;
        } else {
            boolean             hasA = tm.indexOf('A') >= 0;
            boolean             hasQ = tm.indexOf('Q') >= 0;
            boolean             hasZ = tm.indexOf('Z') >= 0;
            if ( hasA ) {
                return false;
            }
            else if ( hasZ && (fillTrimmedReads || fillTrimmedReadsZ) ) {
                return true;
            }
            else if ( hasQ && (fillTrimmedReads || fillTrimmedReadsQ) ) {
                return true;
            } else {
                return false;
            }
        }
    }

    private void buildExtendedRef(final ScoredHaplotype scoredHaplotype, final ReferenceDataSource ref, final SimpleInterval loc, final GATKRead read) {

        // assume no extension
        int     extendStart = 0;
        int     extendEnd = 0;

        // calc soft extension
        if ( fillSoftclippedReads ) {
            final CigarElement elem = !read.isReverseStrand()
                    ? read.getCigar().getLastCigarElement() : read.getCigar().getFirstCigarElement();
            if (elem.getOperator() == CigarOperator.S) {
                if (!read.isReverseStrand()) {
                    extendEnd += elem.getLength();
                } else {
                    extendStart += elem.getLength();
                }
            }
        }

        // add padding
        if ( !read.isReverseStrand() ) {
            extendEnd += haplotypeOutputPaddingSize;
        } else {
            extendStart += haplotypeOutputPaddingSize;
        }

        // add extra fill from haplotype
        if ( (outputFlowLength != 0) && shouldFillFromHaplotype(read) ) {
            int     length = (loc.getEnd() + extendEnd) - (loc.getStart() - extendStart);
            int     delta = Math.max(0, outputFlowLength - length) + EXTRA_FILL_FROM_HAPLOTYPE;
            if ( !read.isReverseStrand() ) {
                extendEnd += delta;
            } else {
                extendStart += delta;
            }
        }

        // compensate for skipage of first hmer
        final int delta = scoredHaplotype.ref.getInterval().size() - scoredHaplotype.haplotype.getGenomeLocation().getLengthOnReference();
        if ( delta != 0 ) {
            if ( !read.isReverseStrand() ) {
                extendStart -= delta;
            } else {
                extendEnd -= delta;
            }
        }


        int         minStart = ref.getSequenceDictionary().getSequence(loc.getContig()).getStart();
        int         maxEnd = ref.getSequenceDictionary().getSequence(loc.getContig()).getEnd();
        if ( isStartSoftclipped(read) ) {
            scoredHaplotype.clippedRef = new ReferenceContext(ref,
                    new SimpleInterval(loc.getContig(), Math.max(minStart, loc.getStart() - extendStart), Math.min(maxEnd, loc.getEnd() + extendEnd)));
        }

        // add front unclipped
        {
            final CigarElement frontElem = !read.isReverseStrand()
                    ? read.getCigar().getFirstCigarElement() : read.getCigar().getLastCigarElement();
            if (frontElem.getOperator() == CigarOperator.S) {
                if (!read.isReverseStrand()) {
                    extendStart += frontElem.getLength();
                } else {
                    extendEnd += frontElem.getLength();
                }
            }
            scoredHaplotype.unclippedRef = new ReferenceContext(ref,
                    new SimpleInterval(loc.getContig(), Math.max(minStart, loc.getStart() - extendStart), Math.min(maxEnd, loc.getEnd() + extendEnd)));
        }
    }

    private boolean arsSame(final Haplotype h1, final Haplotype h2, boolean isReversed) {

        return Arrays.equals(h1.getBases(), h2.getBases());
    }

    private boolean areSame(final Haplotype h, final ReferenceContext ref, boolean isReversed) {
        return Arrays.equals(h.getBases(), reverseComplement(ref.getBases(), isReversed));
    }

    private FlowBasedRead buildFlowRead(final GATKRead read) {

        ReadGroupInfo rgInfo = getReadGroupInfo(getHeaderForReads(), read);

        return new FlowBasedRead(read, rgInfo.flowOrder, rgInfo.maxClass, hcArgs.fbargs);
    }

    private boolean isEndSoftclipped(final GATKRead read) {

        if ( !read.isReverseStrand() ) {
            return read.getCigar().getLastCigarElement().getOperator() == CigarOperator.S;
        } else {
            return read.getCigar().getFirstCigarElement().getOperator() == CigarOperator.S;
        }
    }

    private boolean isStartSoftclipped(final GATKRead read) {

        if ( !read.isReverseStrand() ) {
            return read.getCigar().getFirstCigarElement().getOperator() == CigarOperator.S;
        } else {
            return read.getCigar().getLastCigarElement().getOperator() == CigarOperator.S;
        }
    }

    private boolean isEndPolyTSoftclipped(final GATKRead read) {

        // must be softclipped
        if ( !isEndSoftclipped(read) )
            return false;

        // are all softclipped bases T
        final byte[]      bases = read.getBasesNoCopy();
        if ( !read.isReverseStrand() ) {
            final int      length = read.getCigar().getFirstCigarElement().getLength();
            for ( int n = 0 ; n < length ; n++ ) {
                if (bases[n] != 'T') {
                    return false;
                }
            }
        } else {
            final int      length = read.getCigar().getLastCigarElement().getLength();
            for ( int n = 0 ; n < length ; n++ ) {
                if (bases[bases.length - n - 1] != 'A') {
                    return false;
                }
            }
        }

        // if here all softclipped bases are T
        return true;
    }

    private void debugLog(final GATKRead read, final ReferenceContext referenceContext, final ScoredHaplotype maternal, final ScoredHaplotype paternal) {

        if ( debugMode ) {
            logger.info("read: " + read.getName() + " " + read.getCigar() + " " + read.getFlags());
            logger.info("read:          " + new SimpleInterval(read) + " " + new String(read.getBases()));
            logger.info("ref:           " + new SimpleInterval(referenceContext) + " " + new String(referenceContext.getBases()));
            logger.info("mRef: " + maternal.ref.getInterval() + " " + new String(maternal.ref.getBases()));
            logger.info("pRef: " + paternal.ref.getInterval() + " " + new String(paternal.ref.getBases()));
            logger.info("pmDiff:                                 " + new String(debugBinDiff(maternal.ref.getBases(), paternal.ref.getBases())));
            logger.info("mHap: " + maternal.score + " " + maternal.haplotype);
            logger.info("pHap: " + paternal.score + " " + paternal.haplotype);
        }
    }

    private byte[] debugBinDiff(final byte[] b1, final byte[] b2) {
        final int     len = Math.min(b1.length, b2.length);
        final byte[]  result = new byte[len];

        for ( int n = 0 ; n < len ; n++ ) {
            result[n] = (b1[n] == b2[n]) ? (byte)'_' : (byte)'1';
        }

        return result;

    }

    private double getFlowBasedReadQuality(final FlowBasedRead read, final int maxClass) {

        double      sum = 0;
        for ( int n = 0 ; n < read.getKeyLength() ; n++ ) {
            sum += read.getProb(n, maxClass);
        }
        return sum;
    }

    private Haplotype buildReferenceHaplotype(final ReferenceContext ref, final GATKRead read) {

        Locatable loc = new SimpleInterval(ref.getInterval());
        byte[]          haplotypeBases = reverseComplement(ref.getBases(), read.isReverseStrand());
        final byte[]    readBases = reverseComplement(getSoftclippedBases(read), read.isReverseStrand());

        // skip haplotype bases until same base as read starts (false SNP compensation)
        if ( haplotypeBases[0] != readBases[0] ) {
            final int       skip = detectFalseSNP(haplotypeBases, readBases);
            if ( skip != 0 ) {
                haplotypeBases = Arrays.copyOfRange(haplotypeBases, skip, haplotypeBases.length);
                if ( !read.isReverseStrand() ) {
                    loc = new SimpleInterval(loc.getContig(), loc.getStart() + skip, loc.getEnd());
                } else {
                    loc = new SimpleInterval(loc.getContig(), loc.getStart(), loc.getEnd() - skip);
                }
            }

        }

        final Haplotype haplotype = new Haplotype(haplotypeBases, loc);
        haplotype.setCigar(new CigarBuilder(false)
                .add(new CigarElement(haplotype.length(), CigarOperator.M)).make());

        return haplotype;
    }

    private byte[] getSoftclippedBases(final GATKRead read) {

        final int   startClip = (read.getCigar().getFirstCigarElement().getOperator() == CigarOperator.SOFT_CLIP)
                                    ? read.getCigar().getFirstCigarElement().getLength() : 0;
        final int   endClip = (read.getCigar().getLastCigarElement().getOperator() == CigarOperator.SOFT_CLIP)
                ? read.getCigar().getLastCigarElement().getLength() : 0;

        final byte[] bases = read.getBasesNoCopy();
        if ( startClip == 0 && endClip == 0 ) {
            return bases;
        } else {
            return Arrays.copyOfRange(bases, startClip, bases.length - endClip);
        }
    }

    private int detectFalseSNP(final byte[] haplotypeBases, final byte[] readBases) {

        // parametrisation
        final int       maxSkipageSize = 5;     // will not skip more than this
        final int       requiredRemainingMatchSize = 5; // after skipage, this number of bases must match

        // this might be redundant, but just in case
        if ( haplotypeBases[0] == readBases[0] )
            return 0;

        // establish sizes of homopolymer on read
        int             readHomoSize = 0;
        for ( ; readHomoSize < readBases.length ; readHomoSize++ )
            if ( readBases[readHomoSize] != readBases[0] )
                break;

        // loop on possible skipage
        for ( int skip = 1 ; skip <= maxSkipageSize ; skip++ ) {

            // there must be enough bases
            if ( skip + requiredRemainingMatchSize > haplotypeBases.length ) {
                break;
            }
            if ( Math.max(skip + 1, requiredRemainingMatchSize) > readBases.length ) {
                break;
            }
            // skipage + 1 must be inside homo polymer
            if ( skip + 1 > readHomoSize ) {
                break;
            }

            // remaining bases must match
            if ( arrayRangeEquals(readBases, skip, requiredRemainingMatchSize, haplotypeBases, skip) ) {
                // found skipape point
                return skip;
            }
        }

        // if here, did not find skipage point
        return 0;
    }

    private boolean arrayRangeEquals(final byte[] a1, final int ofs1, final int len, final byte[] a2, final int ofs2) {

        for ( int i = 0 ; i < len ; i++ ) {
            if (a1[ofs1 + i] != a2[ofs2 + i]) {
                return false;
            }
        }
        return true;
    }

    private double scoreReadAgainstHaplotype(final GATKRead read, final ScoredHaplotype sh) {

        // build haplotypes
        final ReadGroupInfo rgInfo = getReadGroupInfo(getHeaderForReads(), read);
        final FlowBasedHaplotype flowHaplotype = new FlowBasedHaplotype(sh.haplotype, rgInfo.flowOrder);

        // create flow read
        final FlowBasedRead   flowRead = new FlowBasedRead(read, rgInfo.flowOrder, rgInfo.maxClass, hcArgs.fbargs);
        if ( read.isReverseStrand() ) {
            flowRead.setDirection(FlowBasedRead.Direction.SYNTHESIS);
            flowRead.applyAlignment();
        }

        if ( !flowRead.isValid() ) {
            return -1;
        }

        // compute alternative score
        final int         hapKeyLength = flowHaplotype.getKeyLength();
        final double      score = FlowFeatureMapper.computeLikelihoodLocal(flowRead, flowHaplotype, hapKeyLength, false);

        return score;
    }

    private double scoreReadAgainstReference(final GATKRead read, final ReferenceContext ref) {

        // build haplotypes
        final ReadGroupInfo           rgInfo = getReadGroupInfo(getHeaderForReads(), read);
        final FlowBasedHaplotype      flowHaplotype = new FlowBasedHaplotype(buildReferenceHaplotype(ref, read), rgInfo.flowOrder);

        // create flow read
        final FlowBasedRead           flowRead = new FlowBasedRead(read, rgInfo.flowOrder, rgInfo.maxClass, hcArgs.fbargs);
        if ( read.isReverseStrand() ) {
            flowRead.setDirection(FlowBasedRead.Direction.SYNTHESIS);
            flowRead.applyAlignment();
        }

        if ( !flowRead.isValid() ) {
            return -1;
        }

        // compute alternative score
        final int         hapKeyLength = flowHaplotype.getKeyLength();
        final double      score = FlowFeatureMapper.computeLikelihoodLocal(flowRead, flowHaplotype, hapKeyLength, false);

        // debug
        if ( debugMode ) {
            logger.info("flowRead: " + flowRead);
            logger.info("flowHaplotype: " + flowHaplotype);
            logger.info("flowRead.key:      " + Arrays.toString(flowRead.getKey()));
            logger.info("flowHaplotype.key: " + Arrays.toString(flowHaplotype.getKey()));
            logger.info("scoreReadAgainstReference: score: " + score);
        }

        return score;
    }

    private int[] buildHaplotypeKey(final String haplotypeSeq, final ReadGroupInfo rgInfo, final boolean isReversed) {

        // create a haplotype to contain the sequence
        final byte[]              seq = reverseComplement(haplotypeSeq.getBytes(), isReversed);
        final Haplotype           h = new Haplotype(seq);
        final FlowBasedHaplotype  flowHaplotype = new FlowBasedHaplotype(h, !isReversed ? rgInfo.flowOrder : rgInfo.reversedFlowOrder);

        // need to start on a T - find out T offset on the flow order
        int[]                     hapKey = flowHaplotype.getKey();
        byte[]                    hapFlowOrder = flowHaplotype.getFlowOrderArray();
        int                       appendZeroCount = 0;
        while ( hapKey[0] == 0 ) {
            hapKey = Arrays.copyOfRange(hapKey, 1, hapKey.length);
        }
        if ( (seq[0] != 'T') && (seq[0] != 'N') ) {
            int     ofs = 0;
            while ( hapFlowOrder[ofs] != 'T' )
                ofs++;
            while ( hapFlowOrder[ofs] != seq[0] ) {
                appendZeroCount++;
                ofs = (ofs + 1) % hapFlowOrder.length;
            }
        }

        if ( appendZeroCount == 0 ) {
            return hapKey;
        } else {
            int[]       hapKey1 = new int[appendZeroCount + hapKey.length];
            System.arraycopy(hapKey, 0, hapKey1, appendZeroCount, hapKey.length);
            return hapKey1;
        }

    }

    private int[] buildHaplotypeKeyForOutput(ScoredHaplotype scoredHaplotype, final ReadGroupInfo rgInfo, final int fillValue, final GATKRead read) {

        boolean                   isReversed = read.isReverseStrand();

        // create key from filled and unclipped version
        int[]                     hapKey = buildHaplotypeKey(new String(scoredHaplotype.unclippedRef.getBases()), rgInfo, isReversed);
        if ( isStartSoftclipped(read) ) {
            int[] hapKeyClipped = buildHaplotypeKey(new String(scoredHaplotype.clippedRef.getBases()), rgInfo, isReversed);
            scoredHaplotype.softclipFrontFillCount = hapKey.length - hapKeyClipped.length;
        } else {
            scoredHaplotype.softclipFrontFillCount = 0;
        }

        // prepare key
        final int flowLength = (outputFlowLength != 0) ? outputFlowLength : hapKey.length;
        final int[] key = new int[flowLength];
        int ofs;
        System.arraycopy(hapKey, 0, key, 0, ofs = Math.min(flowLength, hapKey.length));

        // adjust to a fixed length
        for (  ; ofs < flowLength ; ofs++ ) {
            key[ofs] = fillValue;
        }

        return key;
    }

    private String buildHaplotypeSequenceForOutput(final ScoredHaplotype haplotype, final boolean isReversed, final int keyBaseCount) {

        final StringBuilder      sb = new StringBuilder();
        if ( prependSequence != null ) {
            sb.append(prependSequence);
        }

        final String            seq = new String(reverseComplement(haplotype.unclippedRef.getBases(), isReversed));
        final String            baseCountSeq = seq.substring(0, keyBaseCount);
        sb.append(baseCountSeq);

        if ( appendSequence != null ) {
            sb.append(appendSequence);
        }

        return sb.toString();
    }

    private int[] buildConsensusKey(final int[] k1, final int[] k2) {

        final int         len = Math.min(k1.length, k2.length);
        final int[]       key = new int[len];

        for ( int n = 0 ; n < len ; n++ ) {
            key[n] = (k1[n] == k2[n]) ? k1[n] : -72;
        }

        return key;
    }


    private String flowKeyAsCsvString(final int[] key) {
        return "\"" + Arrays.toString(key).replaceAll("\\[|\\]|\\s", "") + "\"";
    }

    private String flowKeyAsCsvString(int[] key, final String seq, final String flowOrder) {
        final StringBuilder     sb = new StringBuilder();

        sb.append("\"");

        while ( key[0] == 0 ) {
            key = Arrays.copyOfRange(key, 1, key.length);
        }
        if ( (seq.charAt(0) != 'T') && (seq.charAt(0) != 'N') ) {
            int     ofs = 0;
            while ( flowOrder.charAt(ofs) != 'T' )
                ofs++;
            while ( flowOrder.charAt(ofs) != seq.charAt(0) ) {
                sb.append("0,");
                ofs = (ofs + 1) % flowOrder.length();
            }
        }

        sb.append(Arrays.toString(key).replaceAll("\\[|\\]|\\s", ""));

        sb.append("\"");

        return sb.toString();
    }

    private void emitCsvHeaders() {

        outputCsv.println(StringUtils.join(CSV_FIELD_ORDER, ","));
    }

    private void emit(final GATKRead read, final FlowBasedRead flowRead, final double refScore, final ScoredHaplotype maternal, final ScoredHaplotype paternal) throws IOException {

        // build line columns
        final Map<String,Object>        cols = new LinkedHashMap<>();

        // establish fill value
        final String        tm = read.getAttributeAsString("tm");
        boolean             hasA = (tm != null) && tm.indexOf('A') >= 0;
        boolean             hasQ = (tm != null) && tm.indexOf('Q') >= 0;
        boolean             hasZ = (tm != null) && tm.indexOf('Z') >= 0;
        int                 fillValue;
        if ( isEndSoftclipped(read) )
            fillValue = SOFTCLIP_FILL_VALUE;
        else if ( hasQ || hasZ ) {
            fillValue = hasA ? UNKNOWN_FILL_VALUE : NONREF_FILL_VALUE;
        } else {
            fillValue = DEFAULT_FILL_VALUE;
        }

        // read name
        cols.put("ReadName", read.getName());

        // haplotypes and reference scores
        cols.put("PaternalHaplotypeScore", paternal.score);
        cols.put("MaternalHaplotypeScore", maternal.score);
        cols.put("RefHaplotypeScore", refScore);

        // build haplotype keys
        final ReadGroupInfo rgInfo = getReadGroupInfo(getHeaderForReads(), read);
        final int[]           paternalHaplotypeKey = buildHaplotypeKeyForOutput(paternal, rgInfo,fillValue, read);
        final int[]           maternalHaplotypeKey = buildHaplotypeKeyForOutput(maternal, rgInfo,fillValue, read);

        // build haplotype sequence
        final String           paternalHaplotypeSeq = buildHaplotypeSequenceForOutput(paternal, read.isReverseStrand(), keyBases(paternalHaplotypeKey));
        final String           maternalHaplotypeSeq = buildHaplotypeSequenceForOutput(maternal, read.isReverseStrand(), keyBases(maternalHaplotypeKey));

        // fill softclip at front
        softclipFill(paternal, paternalHaplotypeKey);
        softclipFill(maternal, maternalHaplotypeKey);

        // select best and establish consensus
        final boolean          ancestralHaplotypesSame = paternalHaplotypeSeq.equals(maternalHaplotypeSeq);
        final ScoredHaplotype  bestHaplotype = (paternal.score > maternal.score) ? paternal: maternal;
        final int[]            bestHaplotypeKey = (bestHaplotype == paternal) ? paternalHaplotypeKey : maternalHaplotypeKey;
        final int[]            consensus = buildConsensusKey(paternalHaplotypeKey, maternalHaplotypeKey);

        // emit best haplotype
        cols.put("BestHaplotypeSequence", (bestHaplotype == paternal) ? paternalHaplotypeSeq : maternalHaplotypeSeq);
        if ( !ancestralHaplotypesSame )
            cols.put("BestHaplotypeKey", flowKeyAsCsvString(bestHaplotypeKey));
        else
            cols.put("BestHaplotypeKey", flowKeyAsCsvString(consensus));

        // write consensus haplotype
        cols.put("ConsensusHaplotypeKey", flowKeyAsCsvString(consensus));

        // additional  fields
        cols.put("ReadChrom", read.getContig());
        cols.put("ReadStart", read.getStart());
        cols.put("ReadEnd", read.getEnd());
        cols.put("ReadUnclippedStart", read.getUnclippedStart());
        cols.put("ReadUnclippedEnd", read.getUnclippedEnd());
        cols.put("ReadCigar", read.getCigar());

        final String    readSeq = reverseComplement(read.getBasesString(), read.isReverseStrand());
        final int[]     readKey = reverse(flowRead.getKey(), read.isReverseStrand());
        final String    readFlowOrder = reverseComplement(getReadGroupInfo(getHeaderForReads(), read).flowOrder, read.isReverseStrand());
        cols.put("ReadSequence", readSeq);
        cols.put("ReadKey", flowKeyAsCsvString(readKey, readSeq, readFlowOrder));
        cols.put("PaternalHaplotypeInterval", paternal.ref.getInterval());
        cols.put("PaternalHaplotypeSequence", paternalHaplotypeSeq);
        cols.put("MaternalHaplotypeInterval", maternal.ref.getInterval());
        cols.put("MaternalHaplotypeSequence", maternalHaplotypeSeq);

        cols.put("tm", (read.hasAttribute("tm") ? read.getAttributeAsString("tm") : ""));
        cols.put("mapq", read.getMappingQuality());
        cols.put("flags", read.getFlags());

        // construct line
        StringBuilder       sb = new StringBuilder();
        int                 colIndex = 0;
        for ( String field : CSV_FIELD_ORDER ) {
            if ( colIndex++ > 0 ) {
                sb.append(',');
            }
            if ( !cols.containsKey(field) ) {
                throw new GATKException("column missing from csv line: " + field);
            }
            sb.append(cols.get(field));
            cols.remove(field);
        }
        if ( cols.size() > 0 ) {
            throw new GATKException("invalid columns on csv line: " + cols.keySet());
        }

        // output line
        if ( !noOutput ) {
            outputCsv.println(sb);
        }
    }

    private void softclipFill(ScoredHaplotype scoredHaplotype, int[] key) {
        if ( !fillSoftclippedReads ) {
            int     limit = Math.min(scoredHaplotype.softclipFrontFillCount, key.length);
            for (int n = 0; n < limit ; n++) {
                key[n] = SOFTCLIP_FILL_VALUE;
            }
        }
    }

    private int keyBases(int[] key) {
        int     count = 0;
        for ( int c : key ) {
            if (c > 0) {
                count += c;
            }
        }
        return count;
    }

    private synchronized ReadGroupInfo getReadGroupInfo(SAMFileHeader headerForReads, GATKRead read) {

        String              rg = read.getReadGroup();
        ReadGroupInfo       info = readGroupInfo.get(rg);
        if ( info == null ) {
            info = new ReadGroupInfo();
            info.flowOrder = headerForReads.getReadGroup(rg).getFlowOrder();
            info.reversedFlowOrder = reverseComplement(info.flowOrder);
            String mc_string = headerForReads.getReadGroup(rg).getAttribute("mc");
            if ( mc_string == null ) {
                info.maxClass = FlowBasedRead.MAX_CLASS;
            } else {
                info.maxClass = Integer.parseInt(mc_string);
            }
            readGroupInfo.put(rg, info);
        }
        return info;
    }

    private byte[] reverseComplement(final byte[] bases) {

        final byte[] result = new byte[bases.length];
        System.arraycopy(bases, 0, result, 0, result.length);
        SequenceUtil.reverseComplement(result);

        return result;
    }

    private byte[] reverseComplement(final byte[] bases, final boolean isReversed) {
        return !isReversed ? bases : reverseComplement(bases);
    }

    private String reverseComplement(final String bases) {
        return new String(reverseComplement(bases.getBytes()));
    }

    private String reverseComplement(final String bases, final boolean isReversed) {
        return !isReversed ? bases : reverseComplement(bases);
    }

    private int[] reverse(final int[] bytes) {
        final int[] result = new int[bytes.length];
        System.arraycopy(bytes, 0, result, 0, result.length);
        FlowBasedRead.reverse(result, result.length);

        return result;
    }

    private int[] reverse(final int[] bytes, final boolean isReversed) {
        return !isReversed ? bytes : reverse(bytes);
    }

}

