package org.broadinstitute.hellbender.engine.filters.flow;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class HmerQualitySymetricReadFilter extends FlowBasedTPAttributeSymetricReadFilter {
    private static final long serialVersionUID = 1l;

    public HmerQualitySymetricReadFilter() {

    }

    public HmerQualitySymetricReadFilter(final SAMFileHeader header) {
        setHeader(header);
    }

    @Override
    protected byte[] getValuesOfInterest(final GATKRead read) {
        return read.getBaseQualitiesNoCopy();
    }


}
