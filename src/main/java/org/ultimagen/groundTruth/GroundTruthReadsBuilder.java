package org.ultimagen.groundTruth;

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
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.ultimagen.featureMapping.FlowFeatureMapper;
import org.ultimagen.flowBasedRead.alignment.FlowBasedAlignmentEngine;
import org.ultimagen.flowBasedRead.read.FlowBasedHaplotype;
import org.ultimagen.flowBasedRead.read.FlowBasedRead;
import picard.cmdline.programgroups.BaseCallingProgramGroup;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Random;

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

    // locals
    final Random                        random = new Random();
    int                                 outputReadsCount = 0;
    ReferenceDataSource                 maternalReference;
    ReferenceDataSource                 paternalReference;
    AncestralContigLocationTranslator   locationTranslator;
    FlowBasedAlignmentEngine            likelihoodCalculationEngine;
    PrintWriter                         outputCsv;
    private final Map<String, ReadGroupInfo> readGroupInfo = new LinkedHashMap<>();


    private static class ScoredHaplotype {
        ReferenceContext        ref;
        ReferenceContext        extendedRef;
        ReferenceContext        filledRef;
        ReferenceContext        unclippedRef;
        int                     softclipFrontFillCount;
        Haplotype               haplotype;
        FlowBasedRead           flowRead;
        FlowBasedHaplotype      flowHaplotype;
        double                  score;

        ScoredHaplotype(FlowBasedRead flowRead) {
            this.flowRead = flowRead;
        }
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
        ReadLikelihoodCalculationEngine engine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(hcArgs.likelihoodArgs, hcArgs.fbargs, false);
        if ( engine instanceof FlowBasedAlignmentEngine) {
            likelihoodCalculationEngine = (FlowBasedAlignmentEngine)engine;
        } else {
            throw new GATKException("must use a flow based likelihood calculation engine");
        }

        // open output, write header
        outputCsv = new PrintWriter(outputCsvPath.getOutputStream());
        emitCsvHeaders();
    }

    @Override
    public void closeTool() {

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
        if ( discardNonPolytSoftclippedReads && isSoftclipped(read) && !isPolyTSoftclipped(read) ) {
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
            final ScoredHaplotype                   maternal = new ScoredHaplotype(flowRead);
            final ScoredHaplotype                   paternal = new ScoredHaplotype(flowRead);

            // translate location to ascentors
            final Tuple<SimpleInterval, SimpleInterval>  ancestralLocs = locationTranslator.translate(read);
            maternal.ref = new ReferenceContext(maternalReference, ancestralLocs.a);
            paternal.ref = new ReferenceContext(paternalReference, ancestralLocs.b);
            buildExtendedRef(maternal, maternalReference, ancestralLocs.a, read);
            buildExtendedRef(paternal, paternalReference, ancestralLocs.b, read);

            // build haplotypes
            maternal.haplotype = buildReferenceHaplotype(maternal.ref);
            paternal.haplotype = buildReferenceHaplotype(paternal.ref);

            // generate score for reference
            final double        refScore = scoreReadAgainstReference(read, referenceContext);

            // score read against haplotypes, create flow versions of read nad haplotype
            if ( areSame(maternal.haplotype, referenceContext) ) {
                maternal.score = refScore;
            } else {
                maternal.score = scoreReadAgainstHaplotype(read, maternal);
            }
            if ( areSame(paternal.haplotype, referenceContext) ) {
                paternal.score = refScore;
            } else if ( arsSame(maternal.haplotype, paternal.haplotype) ) {
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

        } catch (IOException e) {
            throw new GATKException("failed to process read: " + read.getName(), e);
        }
    }

    private boolean shouldFillFromHaplotype(final GATKRead read) {

        // softclip has priori
        if ( isSoftclipped(read) )
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

        // extend location and build
        scoredHaplotype.extendedRef = new ReferenceContext(ref,
                new SimpleInterval(loc.getContig(), loc.getStart() - extendStart, loc.getEnd() + extendEnd));

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
        scoredHaplotype.filledRef = new ReferenceContext(ref,
                new SimpleInterval(loc.getContig(), loc.getStart() - extendStart, loc.getEnd() + extendEnd));

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
                    new SimpleInterval(loc.getContig(), loc.getStart() - extendStart, loc.getEnd() + extendEnd));
        }
    }

    private boolean arsSame(final Haplotype h1, final Haplotype h2) {

        return Arrays.equals(h1.getBases(), h2.getBases());
    }

    private boolean areSame(final Haplotype h, final ReferenceContext ref) {
        return Arrays.equals(h.getBases(), ref.getBases());
    }

    private FlowBasedRead buildFlowRead(final GATKRead read) {

        ReadGroupInfo rgInfo = getReadGroupInfo(getHeaderForReads(), read);

        return new FlowBasedRead(read, rgInfo.flowOrder, rgInfo.maxClass, hcArgs.fbargs);
    }

    private boolean isSoftclipped(final GATKRead read) {

        if ( !read.isReverseStrand() ) {
            return read.getCigar().getLastCigarElement().getOperator() == CigarOperator.S;
        } else {
            return read.getCigar().getFirstCigarElement().getOperator() == CigarOperator.S;
        }
    }

    private boolean isPolyTSoftclipped(final GATKRead read) {

        // must be softclipped
        if ( !isSoftclipped(read) )
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

    private String getReadNumber(final GATKRead read) {

        final String      name = read.getName();
        final int         i = name.lastIndexOf('-');

        return (i < 0) ? name : name.substring(i + 1);
    }

    private Haplotype buildReferenceHaplotype(final ReferenceContext ref) {

        final Locatable loc = new SimpleInterval(ref.getInterval());
        final Haplotype haplotype = new Haplotype(ref.getBases(), loc);
        haplotype.setCigar(new CigarBuilder(false)
                .add(new CigarElement(haplotype.length(), CigarOperator.M)).make());

        return haplotype;
    }

    private double scoreReadAgainstHaplotype(final GATKRead read, final ScoredHaplotype sh) {

        // build haplotypes
        final ReadGroupInfo rgInfo = getReadGroupInfo(getHeaderForReads(), read);
        sh.flowHaplotype = new FlowBasedHaplotype(sh.haplotype, rgInfo.flowOrder);

        // create flow read
        if ( sh.flowRead == null ) {
            sh.flowRead = new FlowBasedRead(read, rgInfo.flowOrder, rgInfo.maxClass, hcArgs.fbargs);
        }

        if ( !sh.flowRead.isValid() ) {
            return -1;
        }

        // compute alternative score
        final int         hapKeyLength = sh.flowHaplotype.getKeyLength();
        final double      score = FlowFeatureMapper.computeLikelihoodLocal(sh.flowRead, sh.flowHaplotype, hapKeyLength, false);

        return score;
    }

    private double scoreReadAgainstReference(final GATKRead read, final ReferenceContext ref) {

        // build haplotypes
        final ReadGroupInfo           rgInfo = getReadGroupInfo(getHeaderForReads(), read);
        final FlowBasedHaplotype      flowHaplotype = new FlowBasedHaplotype(buildReferenceHaplotype(ref), rgInfo.flowOrder);

        // create flow read
        final FlowBasedRead           flowRead = new FlowBasedRead(read, rgInfo.flowOrder, rgInfo.maxClass, hcArgs.fbargs);

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
        if ( seq[0] != 'T' ) {
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

    private int[] buildHaplotypeKeyForOutput(ScoredHaplotype scoredHaplotype, final ReadGroupInfo rgInfo, final int fillValue, final boolean isReversed) {

        // create key from filled and unclipped version
        int[]                     hapKeyClipped = buildHaplotypeKey(new String(scoredHaplotype.filledRef.getBases()), rgInfo, isReversed);
        int[]                     hapKey = buildHaplotypeKey(new String(scoredHaplotype.unclippedRef.getBases()), rgInfo, isReversed);
        scoredHaplotype.softclipFrontFillCount = hapKey.length - hapKeyClipped.length;

        // prepare key
        final int             flowLength = (outputFlowLength != 0) ? outputFlowLength : hapKey.length;
        final int[]           key = new int[flowLength];
        int                   ofs;
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

    private String flowKeyAsCsvString(byte[] key, final String seq, final String flowOrder) {
        final StringBuilder     sb = new StringBuilder();

        sb.append("\"");

        while ( key[0] == 0 ) {
            key = Arrays.copyOfRange(key, 1, key.length);
        }
        if ( seq.charAt(0) != 'T' ) {
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

        final String[]        fields = {
                "ReadNumber",
                "PaternalHaplotypeScore", "MaternalHaplotypeScore",
                "RefHaplotypeScore",
                "BestHaplotypeSeq", "BestHaplotypeKey", "ConsensusHaplotypeKey",

                "ReadChrom",
                "ReadStart", "ReadEnd",
                "ReadUnclippedStart", "ReadUnclippedEnd",
                "ReadCigar", "ReadSeq", "ReadKey",
                "PaternalHaplotypeInterval", "PaternalHaplotypeSequence",
                "MaternalHaplotypeInterval", "MaternalHaplotypeSequence",
                "tm", "mapq", "flags"
        };

        outputCsv.println(StringUtils.join(fields, ","));
    }

    private void emit(final GATKRead read, final FlowBasedRead flowRead, final double refScore, final ScoredHaplotype maternal, final ScoredHaplotype paternal) throws IOException {

        // build line
        final StringBuilder       sb = new StringBuilder();

        // establish fill value
        final String        tm = read.getAttributeAsString("tm");
        boolean             hasA = (tm != null) && tm.indexOf('A') >= 0;
        boolean             hasQ = (tm != null) && tm.indexOf('Q') >= 0;
        boolean             hasZ = (tm != null) && tm.indexOf('Z') >= 0;
        int                 fillValue;
        if ( isSoftclipped(read) )
            fillValue = SOFTCLIP_FILL_VALUE;
        else if ( hasQ || hasZ ) {
            fillValue = hasA ? UNKNOWN_FILL_VALUE : NONREF_FILL_VALUE;
        } else {
            fillValue = DEFAULT_FILL_VALUE;
        }

        // read number
        final String              readNumber = getReadNumber(read);
        sb.append(readNumber);

        // haplotypes and reference scores
        sb.append("," + paternal.score);
        sb.append("," + maternal.score);
        sb.append("," + refScore);

        // build haplotype keys
        final ReadGroupInfo rgInfo = getReadGroupInfo(getHeaderForReads(), read);
        final int[]           paternalHaplotypeKey = buildHaplotypeKeyForOutput(paternal, rgInfo,fillValue, read.isReverseStrand());
        final int[]           maternalHaplotypeKey = buildHaplotypeKeyForOutput(maternal, rgInfo,fillValue, read.isReverseStrand());

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
        sb.append("," + ((bestHaplotype == paternal) ? paternalHaplotypeSeq : maternalHaplotypeSeq));
        if ( !ancestralHaplotypesSame )
            sb.append("," + flowKeyAsCsvString(bestHaplotypeKey));
        else
            sb.append("," + flowKeyAsCsvString(consensus));

        // write consensus haplotype
        sb.append("," + flowKeyAsCsvString(consensus));

        // additional  fields
        sb.append("," + read.getContig());
        sb.append("," + read.getStart());
        sb.append("," + read.getEnd());
        sb.append("," + read.getUnclippedStart());
        sb.append("," + read.getUnclippedEnd());
        sb.append("," + read.getCigar());

        final String    readSeq = reverseComplement(read.getBasesString(), read.isReverseStrand());
        final byte[]    readKey = reverse(flowRead.getKey(), read.isReverseStrand());
        final String    readFlowOrder = reverseComplement(getReadGroupInfo(getHeaderForReads(), read).flowOrder, read.isReverseStrand());
        sb.append("," + readSeq);
        sb.append("," + flowKeyAsCsvString(readKey, readSeq, readFlowOrder));
        sb.append("," + paternal.ref.getInterval());
        sb.append("," + reverseComplement(paternal.haplotype.getBaseString(), read.isReverseStrand()));
        sb.append("," + maternal.ref.getInterval());
        sb.append("," + reverseComplement(maternal.haplotype.getBaseString(), read.isReverseStrand()));

        sb.append("," + (read.hasAttribute("tm") ? read.getAttributeAsString("tm") : ""));
        sb.append("," + read.getMappingQuality());
        sb.append("," + read.getFlags());

        // write
        outputCsv.println(sb);
    }

    private void softclipFill(ScoredHaplotype scoredHaplotype, int[] key) {
        if ( !fillSoftclippedReads ) {
            for (int n = 0; n < scoredHaplotype.softclipFrontFillCount; n++) {
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
                info.maxClass = 12;
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

    private byte[] reverse(final byte[] bytes) {
        final byte[] result = new byte[bytes.length];
        System.arraycopy(bytes, 0, result, 0, result.length);
        FlowBasedRead.reverse(result, result.length);

        return result;
    }

    private byte[] reverse(final byte[] bytes, final boolean isReversed) {
        return !isReversed ? bytes : reverse(bytes);
    }

}

