package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.FlowBasedHaplotype;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.*;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Filters 'informative' reads from the input SAM/BAM/CRAM file to the SAM/BAM/CRAM file.",
        oneLineSummary = "Filters 'informative' reads in the SAM/BAM/CRAM file",
        programGroup = ReadDataManipulationProgramGroup.class
)
@ExperimentalFeature
@WorkflowProperties
public class SelectInformativeReads extends ReadWalker {

    private static final Logger logger = LogManager.getLogger(SelectInformativeReads.class);

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Write output to this file")
    @WorkflowOutput(optionalCompanions={StandardArgumentDefinitions.OUTPUT_INDEX_COMPANION})
    public GATKPath output;
    private SAMFileGATKReadWriter outputWriter;

    // alleles are provided in vcf form
    @Argument(fullName = ReadFilterArgumentDefinitions.ALLELE_FILE_NAME, doc = "vcf file containing alleles")
    public FeatureDataSource<VariantContext> alleles;

    // width of haplotype expansion around the variant location on each side
    @Argument(fullName = "haplotype-expansion-size", doc = "width of haplotype expansion around the variant location on each side", optional = true)
    public int haplotypeExpansionSize = 5;

    // ref/allele to read distance thresholds
    @Argument(fullName = "min-ref-allele-distance", doc = "min ref allele distance")
    Double minRefAlleleDistance;
    @Argument(fullName = "max-abs-ref-score", doc = "max ref score (absolute value)")
    Double maxAbsRefScore;

    @ArgumentCollection
    public FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();

    @Argument(fullName = "debug-reads", doc = "which reads to debug", optional = true)
    List<String> debugReads;

    @Override
    public void onTraversalStart() {
        outputWriter = createSAMWriter(output, true);
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext ) {
        if ( test(read) ) {
            outputWriter.addRead(read);
        }
    }

    @Override
    public void closeTool() {
        if ( outputWriter != null ) {
            outputWriter.close();
        }
    }

    private boolean test(final GATKRead read) {

        boolean debug = isDebugRead(read);
        if ( debug ) { logger.info("processing: " + read.getName() + " " + read.getCigar() + " " + read.getContig() + ":" + read.getStart() + "-" + read.getEnd());}

        // locate variant contexts that fall within this read
        Boolean testResult = null;
        final SimpleInterval interval = new SimpleInterval(read);
        final Iterator<VariantContext> iterator = alleles.query(interval);
        if ( !iterator.hasNext() ) {
            // fail if does not cross any VCs
            //if ( debug ) { logger.info("has no variants");}
            return false;
        }

        // if has VCs, we'll need the read bases/qualities
        final Pair<byte[], byte[]> readData = AlignmentUtils.getBasesAndBaseQualitiesAlignedOneToOne(read);
        final byte[] readBases = readData.getLeft();
        cleanReadBases(readBases);
        //if ( debug ) { logger.info("readBases: " + new String(readBases));}

        // access read group
        final FlowBasedReadUtils.ReadGroupInfo rgInfo = FlowBasedReadUtils.getReadGroupInfo(getHeaderForReads(), read);

        while ( iterator.hasNext() && testResult == null ) {
            final VariantContext vc = iterator.next();
            if ( debug ) { logger.info(" vc: " + vc.getContig() + ":" + vc.getStart() + "-" + vc.getEnd()); }

            // get bytes under the variant context from the read
            final int vcStartOnReadOffset = vc.getStart() - read.getStart();
            if ( vcStartOnReadOffset < 0 ) {
                if ( debug ) { logger.info("  variant starts before read, ignored"); }
                continue;
            }
            final int vcLength = vc.getEnd() - vc.getStart() + 1;
            if ( vcStartOnReadOffset + vcLength > readBases.length ) {
                if ( debug ) { logger.info("  variant ends after read, ignored"); }
                continue;
            }
            final byte[] vcReadBases = Arrays.copyOfRange(readBases, vcStartOnReadOffset, Math.min(readBases.length - 1, vcStartOnReadOffset + vcLength));
            final byte[] vcRefBases =  vc.getReference().getBases();
            if ( debug ) { logger.info("  vcReadBases/vcRefBases: " + new String(vcReadBases) + "/" + new String(vcRefBases)); }

            // if read bases are the same as reference bases than no need to check further
            if ( Arrays.equals(vcReadBases, vcRefBases) ) {
                //if ( debug ) { logger.info("  same as reference"); }
            } else {
                // loop over non-reference alleles
                for (final Allele allele : vc.getAlleles()) {
                    if (allele.isReference()) {
                        continue;
                    }

                    // get allele bases
                    final byte[] alleleBases = allele.getBases();
                    if ( debug ) { logger.info("  checking allele: " + new String(alleleBases)); }

                    // if read data under the allele is same as allele -> pass
                    if (Arrays.equals(vcReadBases, alleleBases)) {
                        if ( debug ) { logger.info("   same as allele"); }
                    } else {

                        // not the same: generate haplotypes around the location
                        final int vcStartOnRead = vc.getStart() - read.getStart();
                        final int vcEndOnRead = vc.getEnd() - read.getStart();

                        // find ends of haplotypes on the read - must be at least N bases and
                        // contain the last HMER in full
                        int leftExpIndex = Math.max(0, vcStartOnRead - haplotypeExpansionSize);
                        while (leftExpIndex > 0 && readBases[leftExpIndex] == readBases[leftExpIndex - 1]) {
                            leftExpIndex--;
                        }
                        int rightExpIndex = Math.min(readBases.length - 1, vcEndOnRead + 1 + haplotypeExpansionSize);
                        while (rightExpIndex < (readBases.length - 1) && readBases[rightExpIndex] == readBases[rightExpIndex + 1]) {
                            rightExpIndex++;
                        }

                        final byte[] prefixBases = Arrays.copyOfRange(readBases, leftExpIndex, vcStartOnRead);
                        final byte[] suffixBases = Arrays.copyOfRange(readBases, Math.min(readBases.length - 1, vcEndOnRead + 1), rightExpIndex);
                        final Haplotype refHaplotpye = makeHaplotype(prefixBases, vc.getReference().getBases(), suffixBases, true, vc.getStart());
                        final Haplotype alleleHaplotpye = makeHaplotype(prefixBases, alleleBases, suffixBases, false, vc.getStart());
                        if ( debug ) { logger.info("   refHaplotpye:    " + refHaplotpye.getBaseString()); }
                        if ( debug ) { logger.info("   alleleHaplotpye: " + alleleHaplotpye.getBaseString()); }

                        // build flow haplotypes
                        final FlowBasedHaplotype refFlowHaplotpye = new FlowBasedHaplotype(refHaplotpye, rgInfo.flowOrder);
                        final FlowBasedHaplotype alleleFlowHaplotpye = new FlowBasedHaplotype(alleleHaplotpye, rgInfo.flowOrder);

                        // create flow read
                        final FlowBasedRead flowRead = new FlowBasedRead(read, rgInfo.flowOrder,
                                rgInfo.maxClass, fbargs);
                        final int diffLeft = leftExpIndex;
                        final int diffRight = readBases.length - rightExpIndex;
                        if ( debug ) { logger.info("   clipLeft/Right: " + diffLeft + "/" + diffRight); }
                        flowRead.applyBaseClipping(Math.max(0, diffLeft), Math.max(diffRight, 0), false);

                        if (!flowRead.isValid()) {
                            if ( debug ) { logger.info("   clipped flow read turned out invalid!"); }
                            return false;
                        }
                        if ( debug ) { logger.info("   clippedRead:     " + debugFlowReadBases(flowRead)); }

                        // compute alternative score
                        final int hapKeyLength = Math.min(refFlowHaplotpye.getKeyLength(), alleleFlowHaplotpye.getKeyLength());
                        final double refScore = FlowFeatureMapper.computeLikelihoodLocal(flowRead, refFlowHaplotpye, hapKeyLength, false);
                        final double alleleScore = FlowFeatureMapper.computeLikelihoodLocal(flowRead, alleleFlowHaplotpye, hapKeyLength, false);
                        if ( debug ) { logger.info("   refScore: " + refScore); }
                        if ( debug ) { logger.info("   alleleScore: " + alleleScore); }

                        // distances must not be too close
                        if (Math.abs(refScore - alleleScore) < minRefAlleleDistance) {
                            if ( debug ) { logger.info("   failing because score are too close"); }
                            testResult = false;
                        }

                        // ref distance must not be too large
                        if (Math.abs(refScore) > maxAbsRefScore) {
                            if ( debug ) { logger.info("   failing because reference score is too weak"); }
                            testResult = false;
                        }
                    }

                    // break out?
                    if (testResult != null) {
                        break;
                    }
                }
            }
        }

        return testResult != null ? testResult : true;
    }

    private String debugFlowReadBases(FlowBasedRead flowRead) {

        final String flowOrder = flowRead.getFlowOrder();
        final StringBuilder sb = new StringBuilder();
        final int[] key = flowRead.getKey();

        int flowIndex = 0;
        for ( int hmer : key ) {
            for ( int i = 0 ; i < hmer ; i++ )
                sb.append(flowOrder.charAt(flowIndex));
            if ( ++flowIndex >= flowOrder.length() )
                flowIndex = 0;
        }

        return sb.toString();
    }

    private boolean isDebugRead(GATKRead read) {
        if ( debugReads == null || debugReads.size() == 0 ) {
            return false;
        } else {
            return debugReads.get(0).equalsIgnoreCase("All") || debugReads.contains(read.getName());
        }
    }

    private void cleanReadBases(final byte[] array) {
        for ( int i = 0 ; i < array.length ; i++ ) {
            switch ( array[i] ) {
                case 'A':
                case 'T':
                case 'C':
                case 'G':
                case 'a':
                case 't':
                case 'c':
                case 'g':
                    break;
                default:
                    array[i] = 'N';
            }
        }
    }

    private Haplotype makeHaplotype(final byte[] prefixBases, final byte[] alleleBases, final byte[] suffixBases,
                                    final boolean isRef, final int start) {

        try {
            ByteArrayOutputStream os = new ByteArrayOutputStream();
            os.write(prefixBases);
            os.write(alleleBases);
            os.write(suffixBases);
            os.close();

            Haplotype hap = new Haplotype(os.toByteArray(), isRef);
            hap.setAlignmentStartHapwrtRef(start);

            return hap;
        } catch (IOException e) {
            throw new GATKException("failed to build haplotype", e);
        }

    }
}
