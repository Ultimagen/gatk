package org.broadinstitute.hellbender.engine.filters.flow;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * A read filter to test if the TP values for each hmer in a flow based read form
 * a polindrome (as they should)
 */
public class FlowBasedTPAttributeSymetricReadFilter extends FlowBasedHmerBasedReadFilter {
    private static final long serialVersionUID = 1l;

    public FlowBasedTPAttributeSymetricReadFilter() {
        super();
    }

    @Override
    protected byte[] getValuesOfInterest(final GATKRead read) {
        return read.getAttributeAsByteArray(FlowBasedRead.FLOW_MATRIX_TAG_NAME);
    }

    @Override
    protected boolean testHmer(final byte[] values, final int hmerStartingOffset, final int hmerLength) {

        // check for symmetry
        return isPalindrome(values, hmerStartingOffset, hmerLength);
    }
}
