package org.broadinstitute.hellbender.engine.filters.flow;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class FlowBasedTPAttributeSymetricReadFilter extends ReadFilter {
    private static final long serialVersionUID = 1l;

    public FlowBasedTPAttributeSymetricReadFilter() {

    }

    public FlowBasedTPAttributeSymetricReadFilter(final SAMFileHeader header) {
        setHeader(header);
    }

    protected byte[] getValuesOfInterest(final GATKRead read) {
        return read.getAttributeAsByteArray("tp");
    }

    protected boolean checkHmer(final byte[] values, final int ofs, final int length) {

        for (int i = 0; i < length; i++) {
            if (values[ofs + i] != values[ofs + length - 1 - i]) {
                return false;
            }
        }

        return true;
    }

    @Override
    public boolean test(final GATKRead read) {

        // access qualities
        final byte[]        values = getValuesOfInterest(read);
        if ( values == null )
            return false;

        // establish if edges are hard clipped
        final boolean       startHardClipped = read.getCigar().getFirstCigarElement().getOperator() == CigarOperator.HARD_CLIP;
        final boolean       endHardClipped = read.getCigar().getLastCigarElement().getOperator() == CigarOperator.HARD_CLIP;

        // iterate over hmers
        final BaseUtils.HmerIterator      iter = new BaseUtils.HmerIterator(read.getBasesNoCopy());
        int     ofs = 0;
        while ( iter.hasNext() ) {

            // find hmer
            final Pair<Byte,Integer>  hmer = iter.next();
            final int                 hmerLength = hmer.getRight();

            // establish first/last
            final boolean             first = ofs == 0;
            final boolean             last = !iter.hasNext();

            // skip edge hmers if hard clipped
            if ( !((first && startHardClipped) || (last && endHardClipped)) ) {
                if (!checkHmer(values, ofs, hmerLength)) {
                    return false;
                }
            }

            // advance
            ofs += hmerLength;
        }

        // if here, all symetric
        return true;
    }
}
