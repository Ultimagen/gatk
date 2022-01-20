package org.broadinstitute.hellbender.engine.filters.flow;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * A read filter to test if the read's readGroup has a flow order associated with it
 */
public class ReadGroupHasFlowOrderReadFilter extends ReadFilter {
    private static final long serialVersionUID = 1l;

    public ReadGroupHasFlowOrderReadFilter() {

    }

    public ReadGroupHasFlowOrderReadFilter(final SAMFileHeader header) {
        setHeader(header);
    }

    @Override
    public boolean test(final GATKRead read) {

        if ( read.getReadGroup() == null )
            return false;
        else if ( samHeader.getReadGroup(read.getReadGroup()) == null )
            return false;
        else if ( samHeader.getReadGroup(read.getReadGroup()).getFlowOrder() == null )
            return false;
        else
            return true;
    }
}
