package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * An implementation of a feature mapper that finds SNPs (SVN)
 *
 * This class only finds SNP that are surrounded by a specific number of bases identical to the reference.
 */

public class SNVMapper extends BaseFeatureMapper implements FeatureMapper {

    public SNVMapper(FlowFeatureMapperArgumentCollection fmArgs) {
        super(fmArgs);
    }

    @Override
    protected FlowFeatureMapper.MappedFeature detectFeature(GATKRead read, ReferenceContext referenceContext, int readOfs, int refOfs) {
        final byte[] bases = read.getBasesNoCopy();
        final byte[] ref = referenceContext.getBases();

        if ( ref[refOfs] != 'N' && (fmArgs.reportAllAlts || (bases[readOfs] != ref[refOfs])) ) {

            // check if surrounded
            boolean surrounded = isSurrounded(read, referenceContext, readOfs, refOfs);
            if ( ignoreBecauseNotSurrounded(surrounded) ) {
                return null;
            }

            // add this feature
            FlowFeatureMapper.MappedFeature feature = FlowFeatureMapper.MappedFeature.makeFeature(
                    FlowFeatureMapperArgumentCollection.MappingFeatureEnum.SNV, read, readOfs, ref[refOfs], referenceContext.getStart() + refOfs, readOfs - refOfs);

            return enrichFeature(feature, surrounded);
        } else {
            return null;
        }
    }

}
