package org.broadinstitute.hellbender.engine.filters.flow;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadTagValueFilter;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class FlowBasedTPAttributeValidReadFilter extends FlowBasedTPAttributeSymetricReadFilter {
    private static final long serialVersionUID = 1l;

    @Argument(fullName = "read-filter-max-hmer",
            doc = "maxHmer to use for testing in the filter", optional = true)
    public int maxHher = 12;

    public FlowBasedTPAttributeValidReadFilter() {

    }

    public FlowBasedTPAttributeValidReadFilter(final SAMFileHeader header) {
        setHeader(header);
    }

    @Override
    protected boolean checkHmer(final byte[] values, final int ofs, final int length) {

        for ( int i = 0 ; i < length ; i++ ) {
            final int      targetValue = values[ofs+i] + length;

            if ( targetValue < 0 || targetValue > maxHher )
                return false;
        }

        return true;
    }
}
