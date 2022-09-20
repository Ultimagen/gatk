package org.broadinstitute.hellbender.engine.filters;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.text.similarity.LevenshteinDistance;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.featuremapping.FlowFeatureMapper;
import org.broadinstitute.hellbender.tools.walkers.featuremapping.FlowFeatureMapperArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.FlowBasedHaplotype;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.FlowBasedReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;

public class AlleleMappingReadFilter extends ReadFilter {

    private static final long serialVersionUID = 1l;

    // alleles are provided in vcf form
    @Argument(fullName = ReadFilterArgumentDefinitions.ALLELE_FILE_NAME, doc = "vcf file containing alleles")
    public FeatureDataSource<VariantContext> alleles;

    // width of haplotype expansion around the variant location on each side
    @Argument(fullName = "haplotype-expansion-size", doc = "width of haplotype expansion around the variant location on each side", optional = true)
    public int haplotypeExpansionSize = 5;

    // ref/allele to read distance thresholds
    @Argument(fullName = "min-ref-allele-distance", doc = "min ref allele distance", optional = true)
    double minRefAlleleDistance = 0.2;
    @Argument(fullName = "max-ref-distance", doc = "max ref distance", optional = true)
    double maxRefDistance = 0.1;

    @ArgumentCollection
    public FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();

    @Override
    public boolean test(GATKRead read) {

        // locate variant contexts that fall within this read
        Boolean testResult = null;
        final SimpleInterval interval = new SimpleInterval(read);
        final Iterator<VariantContext> iterator = alleles.query(interval);
        if ( !iterator.hasNext() ) {
            // fail if does not cross any VCs
            return false;
        }

        // if has VCs, we'll need the read bases/qualities
        final Pair<byte[], byte[]> readData = AlignmentUtils.getBasesAndBaseQualitiesAlignedOneToOne(read);
        final byte[] readBases = readData.getLeft();
        cleanReadBases(readBases);

        // access read group
        final FlowBasedReadUtils.ReadGroupInfo rgInfo = FlowBasedReadUtils.getReadGroupInfo(samHeader, read);

        while ( iterator.hasNext() && testResult == null ) {
            final VariantContext vc = iterator.next();

            // get bytes under the variant context from the read
            final int vcStartOnReadOffset = vc.getStart() - read.getStart();
            if ( vcStartOnReadOffset < 0 ) {
                continue;
            }
            final int vcLength = vc.getEnd() - vc.getStart() + 1;
            if ( vcStartOnReadOffset + vcLength > readBases.length ) {
                continue;
            }
            final byte[] vcReadBases = Arrays.copyOfRange(readBases, vcStartOnReadOffset, Math.min(readBases.length - 1, vcStartOnReadOffset + vcLength));

            // loop over non-reference alleles
            for ( final Allele allele : vc.getAlleles() ) {
                if ( allele.isReference() ) {
                    continue;
                }

                // get allele bases
                final byte[] alleleBases = allele.getBases();

                // if read data under the allele is same as allele -> pass
                if ( Arrays.equals(vcReadBases, alleleBases) ) {
                    ;
                } else {

                    // not the same: generate haplotypes around the location
                    final int vcStartOnRead = vc.getStart() - read.getStart();
                    final int vcEndOnRead = vc.getEnd() - read.getStart();

                    // find ends of haplotypes on the read - must be at least N bases and
                    // contain the last HMER in full
                    // TODO: adjust to HMER boundary
                    int leftExp = haplotypeExpansionSize;
                    int rightExp = haplotypeExpansionSize;

                    final byte[] prefixBases = Arrays.copyOfRange(readBases, Math.max(0, vcStartOnRead - leftExp), vcStartOnRead);
                    final byte[] suffixBases = Arrays.copyOfRange(readBases, Math.min(readBases.length - 1, vcEndOnRead + 1), Math.min(readBases.length - 1, vcEndOnRead + 1 + rightExp));
                    final Haplotype refHaplotpye = makeHaplotype(prefixBases, vc.getReference().getBases(), suffixBases, true, vc.getStart());
                    final Haplotype alleleHaplotpye = makeHaplotype(prefixBases, alleleBases, suffixBases, false, vc.getStart());

                    // build flow haplotypes
                    final FlowBasedHaplotype    refFlowHaplotpye = new FlowBasedHaplotype(refHaplotpye, rgInfo.flowOrder);
                    final FlowBasedHaplotype    alleleFlowHaplotpye = new FlowBasedHaplotype(alleleHaplotpye, rgInfo.flowOrder);

                    // create flow read
                    final FlowBasedRead flowRead = new FlowBasedRead(read, rgInfo.flowOrder,
                            rgInfo.maxClass, fbargs);
                    final int diffLeft = vcStartOnRead;
                    final int diffRight = flowRead.getEnd() - vcEndOnRead;
                    flowRead.applyBaseClipping(Math.max(0, diffLeft), Math.max(diffRight, 0), false);

                    if ( !flowRead.isValid() ) {
                        return false;
                    }

                    // compute alternative score
                    final int         hapKeyLength = Math.min(refFlowHaplotpye.getKeyLength(), alleleFlowHaplotpye.getKeyLength());
                    final double      refScore = FlowFeatureMapper.computeLikelihoodLocal(flowRead, refFlowHaplotpye, hapKeyLength, false);
                    final double      alleleScore = FlowFeatureMapper.computeLikelihoodLocal(flowRead, alleleFlowHaplotpye, hapKeyLength, false);

                    // distances must not be too close
                    if ( Math.abs(refScore - alleleScore) < minRefAlleleDistance ) {
                        testResult = false;
                    }

                    // ref distance must not be too large
                    if ( refScore > maxRefDistance ) {
                        testResult = false;
                    }
                }

                // break out?
                if ( testResult != null ) {
                    break;
                }
            }
        }

        return testResult != null ? testResult : true;
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
