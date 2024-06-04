package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.haplotype.FlowBasedHaplotype;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.FlowBasedReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;

/**
 * An implementation of a feature mapper that finds INDEL features
 *
 * This class only finds features that are surrounded by a specific number of bases identical to the reference.
 */

public class INDELMapper extends BaseFeatureMapper implements FeatureMapper {

    public INDELMapper(FlowFeatureMapperArgumentCollection fmArgs, SAMFileHeader hdr) {
        super(fmArgs, hdr);
    }

    @Override
    protected FlowFeatureMapper.MappedFeature detectFeature(GATKRead read, ReferenceContext referenceContext, final byte bases[], final byte ref[],
                                                            int readOfs, int refOfs, int readLength, int refLength) {

        // check if surrounded
        boolean surrounded = isSurrounded(read, referenceContext, readOfs, refOfs, readLength, refLength);
        if (ignoreBecauseNotSurrounded(surrounded)) {
            return null;
        }

        // check for that it is the kind of INDEL we want
        byte[] readHaplotype = protectedRange(bases, readOfs - surroundBefore, readOfs + readLength + surroundAfter);
        byte[] refHaplotype = protectedRange(ref, refOfs - surroundBefore, refOfs + refLength + surroundAfter);
        if (readHaplotype == null || refHaplotype == null) {
            return null;
        }
        if (!acceptIndel(read, readHaplotype, refHaplotype)) {
            return null;
        }

        // add this feature
        byte[] readBases = protectedRange(bases, readOfs - 1, readOfs + readLength);
        byte[] refBases = protectedRange(ref, refOfs - 1, refOfs + refLength);
        if (readBases == null || refBases == null) {
            return null;
        }

        FlowFeatureMapper.MappedFeature feature = FlowFeatureMapper.MappedFeature.makeFeature(
                FlowFeatureMapperArgumentCollection.MappingFeatureEnum.INDEL,
                read, readBases, refBases, readOfs - 1, referenceContext.getStart() + refOfs - 1, readOfs - refOfs);

        return enrichFeature(feature, surrounded);
    }

    private boolean acceptIndel(GATKRead read, byte[] bases, byte[] ref) {

        // move to flow space
        String flowOrder = new String(FlowBasedReadUtils.getReadFlowOrder(hdr, read));
        FlowBasedHaplotype readHap  = new FlowBasedHaplotype(new Haplotype(bases, false), flowOrder);
        FlowBasedHaplotype refHap = new FlowBasedHaplotype(new Haplotype(ref, false), flowOrder);

        // access flow keys
        int[] readKey = readHap.getKey();
        int[] refKey = refHap.getKey();

        if ( readKey.length != refKey.length ) {
            // this is probably a cycle skip
            // TODO: verify
            return true;
        }

        int changeCount = 0;
        for ( int i = 0 ; i < readKey.length ; i++ ) {
            changeCount += Math.abs(readKey[i] - refKey[i]);
        }

        return changeCount != 1;
    }

    private byte[] protectedRange(final byte[] bases, final int from, final int to) {
        if ( Math.min(from, to) < 0 || Math.max(from, to) > bases.length ) {
            return null;
        } else {
            return Arrays.copyOfRange(bases, from, to);
        }
    }

    @Override
    protected boolean isIndelMapper() {
        return true;
    }

}
