package org.broadinstitute.hellbender.engine.filters;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;
import java.util.Iterator;

public class AlleleMappingReadFilter extends ReadFilter {

    private static final long serialVersionUID = 1l;

    // alleles are provided in vcf form
    @Argument(fullName = ReadFilterArgumentDefinitions.ALLELE_FILE_NAME, doc = "vcf file containing alleles")
    public FeatureDataSource<VariantContext> alleles;

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

        while ( iterator.hasNext() && testResult == null ) {
            final VariantContext vc = iterator.next();

            // get bytes under the variant context from the read
            final int vcStartOnReadOffset = vc.getStart() - read.getStart();
            final int vcLength = vc.getEnd() - vc.getStart() + 1;
            final byte[] vcReadBases = Arrays.copyOfRange(readBases, vcStartOnReadOffset, vcStartOnReadOffset + vcLength);

            // loop over non-reference alleles
            for ( final Allele allele : vc.getAlleles() ) {
                if ( allele.isReference() ) {
                    continue;
                }

                // get allele bases
                final byte[] alleleBases = allele.getBases();

                // NOTE: add actual test here
                if ( Arrays.equals(vcReadBases, alleleBases) ) {
                    testResult = true;
                }

                // break out?
                if ( testResult != null ) {
                    break;
                }
            }
        }

        return testResult != null ? testResult : false;
    }
}
