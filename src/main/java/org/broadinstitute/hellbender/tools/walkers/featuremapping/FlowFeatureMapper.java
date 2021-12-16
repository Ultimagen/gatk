package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.FlowBasedProgramGroup;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriter;
import org.broadinstitute.hellbender.utils.haplotype.FlowBasedHaplotype;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;

import java.util.*;


/**
 * Finds specific features in reads: SNP/indel,
 * scores the confidence of each feature relative to the reference in each read
 * and writes them into a VCF file
 *
 * <p>
 * At this point, this tool finds SNVs
 * </p>
 *
 * <h3> Input </h3>
 * <ul>
 *     <li> Coordinate-sorted and indexed SAM/BAM/CRAM </li>
 * </ul>
 *
 * <h3> Output </h3>
 * <ul>
 *     <li> Coordinate-sorted and indexed VCF </li>
 * </ul>
 *
 * <h3>Usage examples</h3>
 * Find SNVs in chromosome 20.
 * <pre>
 * gatk FlowFeatureMapper \
 *   -I input.bam \
 *   -L 20 \
 *   -O chr20_snv.vcf
 * </pre>
 *
 * {@GATK.walkertype ReadWalker}
 */
@CommandLineProgramProperties(
        summary = "Mapping features (flow space processing)",
        oneLineSummary = "Map/find features in BAM file, output VCF. Initially mapping SNVs",
        programGroup = FlowBasedProgramGroup.class
)


@DocumentedFeature
@ExperimentalFeature
public final class FlowFeatureMapper extends ReadWalker {

    private static final Logger     logger = LogManager.getLogger(FlowFeatureMapper.class);

    private static final String     VCB_SOURCE = "fm";

    private static final String     VCF_READ_NAME = "X_RN";
    private static final String     VCF_SCORE = "X_SCORE";
    private static final String     VCF_FLAGS = "X_FLAGS";
    private static final String     VCF_MAPQ = "X_MAPQ";
    private static final String     VCF_CIGAR = "X_CIGAR";
    private static final String     VCF_READ_COUNT = "X_READ_COUNT";
    private static final String     VCF_FILTERED_COUNT = "X_FILTERED_COUNT";
    private static final String     VCF_FC1 = "X_FC1";
    private static final String     VCF_FC2 = "X_FC2";
    private static final String     VCF_LENGTH = "X_LENGTH";
    private static final String     VCF_EDIST = "X_EDIST";
    private static final String     VCF_INDEX = "X_INDEX";

    private static final Double     LOWEST_PROB = 0.0001;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "File to which variants should be written")
    public GATKPath outputVCF = null;

    @ArgumentCollection
    private FlowFeatureMapperArgumentCollection fmArgs = new FlowFeatureMapperArgumentCollection();

    @ArgumentCollection
    private final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();

    static private class ReadGroupInfo {
        final String  flowOrder;
        final int     maxClass;

        ReadGroupInfo(final String flowOrder, final int maxClass) {
            this.flowOrder = flowOrder;
            this.maxClass = maxClass;
        }
    }

    static private class ReadContext implements Comparable<ReadContext> {
        final GATKRead         read;
        final ReferenceContext referenceContext;

        ReadContext(final GATKRead read, final ReferenceContext referenceContext) {
            this.read = read;
            this.referenceContext = referenceContext;
        }

        @Override
        public int compareTo(ReadContext o) {
            int     delta = read.getContig().compareTo(o.read.getContig());

            delta = (delta != 0) ? delta : Integer.compare(read.getStart(), o.read.getStart());
            delta = (delta != 0) ? delta : Integer.compare(read.getEnd(), o.read.getEnd());

            return delta;
        }
    }

    // locals
    private VariantContextWriter                vcfWriter;
    final private PriorityQueue<Feature>        featureQueue = new PriorityQueue<>();
    final private PriorityQueue<ReadContext>    readQueue = new PriorityQueue<>();
    private FeatureMapper                       mapper;
    private final Map<String, ReadGroupInfo>    readGroupInfo = new LinkedHashMap<>();

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        mapper = buildMapper();

        // open output vcf
        // The HC engine will make the right kind (VCF or GVCF) of writer for us
        final SAMSequenceDictionary sequenceDictionary = getHeaderForReads().getSequenceDictionary();
        vcfWriter = makeVCFWriter(outputVCF, sequenceDictionary, createOutputVariantIndex, createOutputVariantMD5, outputSitesOnlyVCFs);
        vcfWriter.writeHeader(makeVCFHeader(sequenceDictionary, getDefaultToolVCFHeaderLines()));
    }

    @Override
    public void closeTool() {
        flushQueue(null, null);
        super.closeTool();
        vcfWriter.close();
    }

    public VariantContextWriter makeVCFWriter( final GATKPath outputVCF, final SAMSequenceDictionary readsDictionary,
                                               final boolean createOutputVariantIndex, final boolean  createOutputVariantMD5,
                                               final boolean sitesOnlyMode ) {
        Utils.nonNull(outputVCF);
        Utils.nonNull(readsDictionary);

        final List<Options> options = new ArrayList<>(2);
        if (createOutputVariantIndex) {options.add(Options.INDEX_ON_THE_FLY);}
        if (sitesOnlyMode) {options.add(Options.DO_NOT_WRITE_GENOTYPES);}

        VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(
                outputVCF.toPath(),
                readsDictionary,
                createOutputVariantMD5,
                options.toArray(new Options[options.size()])
        );

        if ( hcArgs.emitReferenceConfidence == ReferenceConfidenceMode.GVCF ) {
            try {
                writer = new GVCFWriter(writer, new ArrayList<Number>(hcArgs.GVCFGQBands), hcArgs.floorBlocks);
            } catch ( IllegalArgumentException e ) {
                throw new CommandLineException.BadArgumentValue("GQBands", "are malformed: " + e.getMessage());
            }
        }

        return writer;
    }

    public VCFHeader makeVCFHeader(final SAMSequenceDictionary sequenceDictionary, final Set<VCFHeaderLine>  defaultToolHeaderLines ) {
        final Set<VCFHeaderLine> headerInfo = new HashSet<>();
        headerInfo.addAll(defaultToolHeaderLines);

        // all callers need to add these standard annotation header lines
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerInfo.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));

        // all callers need to add these standard FORMAT field header lines
        VCFStandardHeaderLines.addStandardFormatLines(headerInfo, true,
                VCFConstants.GENOTYPE_KEY,
                VCFConstants.GENOTYPE_QUALITY_KEY,
                VCFConstants.DEPTH_KEY,
                VCFConstants.GENOTYPE_PL_KEY);

        // add our own headers
        headerInfo.add(new VCFInfoHeaderLine(VCF_READ_NAME, 1, VCFHeaderLineType.String, "Read name"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_SCORE, 1, VCFHeaderLineType.Float, "Mapping score"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_FLAGS, 1, VCFHeaderLineType.Integer, "Read flags"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_MAPQ, 1, VCFHeaderLineType.Integer, "Read mapqe"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_CIGAR, 1, VCFHeaderLineType.String, "Read CIGAR"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_READ_COUNT, 1, VCFHeaderLineType.Integer, "Number of reads containing this location"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_FILTERED_COUNT, 1, VCFHeaderLineType.Integer, "Number of reads containing this location that agree with reference according to fitler"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_FC1, 1, VCFHeaderLineType.Integer, "Number of M bases different on read from references"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_FC2, 1, VCFHeaderLineType.Integer, "Number of features before score threshold filter"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_LENGTH, 1, VCFHeaderLineType.Integer, "Read length"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_EDIST, 1, VCFHeaderLineType.Integer, "Read Levenshtein edit distance from reference"));
        headerInfo.add(new VCFInfoHeaderLine(VCF_INDEX, 1, VCFHeaderLineType.Integer, "Ordinal index, from start of the read, where the feature was found"));
        for ( String name : fmArgs.copyAttr ) {
            headerInfo.add(new VCFInfoHeaderLine(fmArgs.copyAttrPrefix + name, 1, VCFHeaderLineType.String, "copy-attr: " + name));
        }

        final VCFHeader vcfHeader = new VCFHeader(headerInfo);
        vcfHeader.setSequenceDictionary(sequenceDictionary);
        return vcfHeader;
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        // include dups?
        if ( read.isDuplicate() && !fmArgs.includeDupReads ) {
            return;
        }

        // include supplementary alignments?
        if ( read.isSupplementaryAlignment() && !fmArgs.keepSupplementaryAlignments ) {
            return;
        }

        // flush qeues up to this read
        flushQueue(read, referenceContext);

        // find features in read
        mapper.forEachOnRead(read, referenceContext, fr -> {
            if ( logger.isDebugEnabled() ) {
                logger.debug("fr: " + fr);
            }

            // score the feature
            fr.score = scoreFeature(fr);

            // emit feature if filters in
            if ( filterFeature(fr) ) {
                featureQueue.add(fr);
            }
        });
    }

    private void flushQueue(final GATKRead read, final ReferenceContext referenceContext) {

        // emit all?
        if ( read == null ) {
            while ( featureQueue.size() != 0 ) {
                final Feature         fr = featureQueue.poll();
                enrichFeature(fr);
                emitFeature(fr);
            }
        } else {
            // enter read into the queue
            readQueue.add(new ReadContext(read, referenceContext));

            // emit all features that start before this read
            while ( featureQueue.size() != 0 ) {
                Feature     fr = featureQueue.peek();
                if ( !fr.read.getContig().equals(read.getContig())
                            || (fr.start < read.getStart()) ) {
                    fr = featureQueue.poll();
                    enrichFeature(fr);
                    emitFeature(fr);
                }
                else {
                    break;
                }
            }

            // remove all reads that start before this read
            while ( readQueue.size() != 0 ) {
                ReadContext rc = readQueue.peek();

                if ( !rc.read.getContig().equals(read.getContig())
                        || (rc.read.getEnd() < read.getStart()) ) {
                    rc = readQueue.poll();
                }
                else {
                    break;
                }
            }
        }
    }

    private void enrichFeature(final Feature fr) {

        // loop on queued reads, count and check if should be counted as filtered
        final Locatable   loc = new SimpleInterval(fr.read.getContig(), fr.start, fr.start);
        for ( ReadContext rc : readQueue ) {
            if ( rc.read.contains(loc) ) {
                fr.readCount++;
                if ( mapper.noFeatureButFilterAt(rc.read, rc.referenceContext, fr.start) ) {
                    fr.filteredCount++;
                }
            }
        }
    }

    private double scoreFeature(final Feature fr) {

        // build haplotypes
        final ReadGroupInfo           rgInfo = getReadGroupInfo(fr.read);
        final FlowBasedHaplotype[]    haplotypes = buildHaplotypes(fr, rgInfo.flowOrder);

        // create flow read
        final FlowBasedRead   flowRead = new FlowBasedRead(fr.read, rgInfo.flowOrder,
                                                                        rgInfo.maxClass, hcArgs.fbargs);
        final int diffLeft = haplotypes[0].getStart() - flowRead.getStart() + fr.offsetDelta;
        final int diffRight = flowRead.getEnd() - haplotypes[0].getEnd();
        flowRead.applyBaseClipping(Math.max(0, diffLeft), Math.max(diffRight, 0), false);

        if ( !flowRead.isValid() ) {
            return -1;
        }

        // compute alternative score
        final int         hapKeyLength = Math.min(haplotypes[0].getKeyLength(), haplotypes[1].getKeyLength());
        final double      readScore = computeLikelihoodLocal(flowRead, haplotypes[0], hapKeyLength, false);
        final double      refScore = computeLikelihoodLocal(flowRead, haplotypes[1], hapKeyLength, false);
        double            score = readScore - refScore;
        if ( !Double.isNaN(fmArgs.limitScore) ) {
            score = Math.min(score, fmArgs.limitScore);
        }

        if ( ((Double.isNaN(score) || (score < 0)) && fmArgs.debugNegatives)
             || (fmArgs.debugReadName != null && fmArgs.debugReadName.contains(fr.read.getName())) ) {
            logger.info("**** debug read: " + fr.read);
            logger.info("readBases: " + fr.read.getBasesString());
            logger.info("flowRead: " + flowRead);
            logger.info("flowBases: " + flowRead.getBasesString());
            logger.info("flowOrder: " + flowRead.getFlowOrder());
            logger.info("flowKey: " + flowRead.getKeyLength() + " " + Arrays.toString(flowRead.getKey()));
            logger.info("readHaplotype: " + haplotypes[0]);
            logger.info("readHapKey: " + haplotypes[0].getKeyLength() + " " + Arrays.toString(haplotypes[0].getKey()));
            computeLikelihoodLocal(flowRead, haplotypes[0], hapKeyLength, true);
            logger.info("refrHaplotype: " + haplotypes[1]);
            logger.info("refrHapKey: " + haplotypes[1].getKeyLength() + " " + Arrays.toString(haplotypes[1].getKey()));
            computeLikelihoodLocal(flowRead, haplotypes[1], hapKeyLength, true);
            logger.info("score: " + score);

            // analyze read
            final FlowBasedRead flowRead2 = new FlowBasedRead(fr.read, rgInfo.flowOrder, rgInfo.maxClass, hcArgs.fbargs);
            final int[]        key2 = flowRead2.getKey();
            for ( int i = 0 ; i < key2.length ; i++ ) {
                final double      p1 = flowRead2.getProb(i, key2[i]);
                for ( int j = 0 ; j < rgInfo.maxClass ; j++ ) {
                    final double      p2 = flowRead2.getProb(i, j);
                    if ( p2 > p1 )
                        logger.info(String.format("prob at %s key[%d]=%d, %f is lower than at %d which is %f",
                                                                flowRead2.getName(), i, key2[i], p1, j, p2));
                }
            }
        }

        if ( score < 0 && !fmArgs.keepNegatives && score != -1.0 ) {
            score = 0;
        }

        return score;
    }

    public static double computeLikelihoodLocal(final FlowBasedRead read, final FlowBasedHaplotype haplotype, final int hapKeyLength, final boolean debug) {

        final byte[] flowOrder = haplotype.getFlowOrderArray();
        final byte   readFlowOrder0 = read.getFlowOrderArray()[0];
        int startingPoint = 0;
        for (int i = 0; i < flowOrder.length; i++) {
            if (flowOrder[i] == readFlowOrder0) {
                startingPoint = i;
                break;
            }
        }
        final int[]         key = haplotype.getKey();

        // debug support
        StringBuffer        debugMessage = null;
        if ( debug )
            debugMessage = new StringBuffer(Integer.toString(startingPoint) + " hmer prob |");
        double              result = 0 ;
        for (int i = 0; i < read.getKeyLength(); i++) {
            int     index = i + startingPoint;
            double  prob = 0;
            int     locationToFetch = 0;
            if ( index < hapKeyLength ) {
                locationToFetch = Math.min(key[index] & 0xff, read.getMaxHmer() + 1);
                prob = read.getProb(i, locationToFetch);
            } else {
                if ( debug ) {
                    debugMessage.append(" clip");
                }
                break;
            }
            if ( prob == 0.0 ) {
                prob = LOWEST_PROB;
            }
            result += Math.log10(prob);

            if ( debug ) {
                debugMessage.append(String.format(" %d %.4f", locationToFetch, prob));
            }
        }

        if ( debug ) {
            debugMessage.append(" | " + result);
            logger.info("debugMessage: " + debugMessage);
        }

        return result;
    }

    private FlowBasedHaplotype[] buildHaplotypes(final Feature fr, final String flowOrder) {

        // build bases for flow haplotypes
        // NOTE!!!: this code assumes length of feature on read and reference is the same
        // this is true for SNP but not for INDELs - it will have to be re-written!
        // TODO: write for INDEL
        final byte[] bases = fr.read.getBasesNoCopy();
        int         offset = fr.readBasesOffset;
        int         refStart = fr.start;
        int         refModOfs = 0;
        if ( offset > 0 ) {
            // reach into hmer before
            offset--;
            refModOfs++;
            refStart--;

            // extend until start of hmer
            final byte        hmerBase = bases[offset];
            while ( offset > 0 && bases[offset-1] == hmerBase ) {
                offset--;
                refModOfs++;
                refStart--;
            }
        }
        final byte[]      sAltBases = Arrays.copyOfRange(bases, offset, bases.length);
        final byte[]      sRefBases = Arrays.copyOf(sAltBases, sAltBases.length);
        System.arraycopy(fr.refBases, 0, sRefBases, refModOfs, fr.refBases.length);

        // construct haplotypes
        final SimpleInterval genomeLoc = new SimpleInterval(fr.read.getContig(), refStart, refStart + sAltBases.length - 1);
        final Cigar          cigar = new Cigar();
        cigar.add(new CigarElement(sAltBases.length, CigarOperator.M));
        final Haplotype      altHaplotype = new Haplotype(sAltBases, false);
        final Haplotype      refHaplotype = new Haplotype(sRefBases, true);
        altHaplotype.setGenomeLocation(genomeLoc);
        refHaplotype.setGenomeLocation(genomeLoc);
        altHaplotype.setCigar(cigar);
        refHaplotype.setCigar(cigar);

        // prepare flow based haplotypes
        final FlowBasedHaplotype[] result = {
                                new FlowBasedHaplotype(altHaplotype, flowOrder),
                                new FlowBasedHaplotype(refHaplotype, flowOrder)
                            };

        // return
        return result;
    }

    private synchronized ReadGroupInfo getReadGroupInfo(final GATKRead read) {

        final String              rg = read.getReadGroup();
        if ( readGroupInfo.containsKey(rg) ) {
            return readGroupInfo.get(rg);
        } else {
            final String mc = getHeaderForReads().getReadGroup(rg).getAttribute("mc");
            final ReadGroupInfo info = new ReadGroupInfo(getHeaderForReads().getReadGroup(rg).getFlowOrder(),
                                            (mc == null) ? FlowBasedRead.MAX_CLASS : Integer.parseInt(mc));
            readGroupInfo.put(rg, info);

            return info;
        }
    }

    private boolean filterFeature(final Feature fr) {

        if ( fmArgs.excludeNaNScores && Double.isNaN(fr.score) ) {
            return false;
        } else if ( fr.score > fmArgs.maxScore ) {
            return false;
        } else if ( fr.score < fmArgs.minScore ) {
            return false;
        }

        return true;
    }

    private void emitFeature(final Feature fr) {

        // create alleles
        final Collection<Allele>          alleles = new LinkedList<>();
        alleles.add(Allele.create(fr.readBases, false));
        alleles.add(Allele.create(fr.refBases, true));

        // create variant context builder
        final VariantContextBuilder       vcb = new VariantContextBuilder(
                                                    VCB_SOURCE,
                                                    fr.read.getContig(),
                                                    fr.start,
                                                    fr.start + fr.refBases.length - 1,
                                                    alleles);

        // copy attributes
        vcb.attribute(VCF_READ_NAME, fr.read.getName());
        vcb.attribute(VCF_SCORE, String.format("%.5f", fr.score));
        vcb.attribute(VCF_FLAGS, fr.read.getFlags());
        vcb.attribute(VCF_MAPQ, fr.read.getMappingQuality());
        vcb.attribute(VCF_CIGAR, fr.read.getCigar().toString());
        vcb.attribute(VCF_READ_COUNT, fr.readCount);
        vcb.attribute(VCF_FILTERED_COUNT, fr.filteredCount);
        vcb.attribute(VCF_FC1, fr.nonIdentMBasesOnRead);
        vcb.attribute(VCF_FC2, fr.featuresOnRead);
        vcb.attribute(VCF_LENGTH, fr.read.getLength());
        vcb.attribute(VCF_EDIST, fr.refEditDistance);
        vcb.attribute(VCF_INDEX, fr.index);
        for ( String name : fmArgs.copyAttr ) {
            if ( fr.read.hasAttribute(name) ) {
                vcb.attribute(fmArgs.copyAttrPrefix + name, fr.read.getAttributeAsString(name));
            }
        }
        final VariantContext      vc = vcb.make();

        // write to file
        vcfWriter.add(vc);
    }

    private FeatureMapper buildMapper() {

        // build appropriate mapper
        if ( fmArgs.mappingFeature == FlowFeatureMapperArgumentCollection.MappingFeatureEnum.SNV ) {
            return new SNVMapper(fmArgs);
        } else {
            throw new GATKException("unsupported mappingFeature: " + fmArgs.mappingFeature);
        }
    }
}

