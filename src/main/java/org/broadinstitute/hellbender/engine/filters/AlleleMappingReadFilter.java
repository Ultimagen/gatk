package org.broadinstitute.hellbender.engine.filters;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Iterator;

public class AlleleMappingReadFilter extends ReadFilter {

    private static final long serialVersionUID = 1l;

    // alleles are provided in vcf form
    @Argument(fullName = ReadFilterArgumentDefinitions.ALLELE_FILE_NAME, doc = "vcf file containing alleles")
    public FeatureDataSource<VariantContext> alleles;

    @Override
    public boolean test(GATKRead read) {

        // locate alleles that fall within this read
        Boolean testResult = null;
        final SimpleInterval interval = new SimpleInterval(read);
        final Iterator<VariantContext> iterator = alleles.query(interval);
        while ( iterator.hasNext() && testResult == null ) {
            final VariantContext vc = iterator.next();

            // for now, we mark the read passed if it crosses any of the alleles
            testResult = true;
        }

        return testResult != null ? testResult : false;
    }
}
