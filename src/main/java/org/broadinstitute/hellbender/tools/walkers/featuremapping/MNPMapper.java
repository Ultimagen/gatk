package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;

/**
 * An implementation of a feature mapper that finds MNP (multi-nucleotide polymorphism) features
 *
 * This class only finds features that are surrounded by a specific number of bases identical to the reference.
 */

public class MNPMapper extends BaseFeatureMapper implements FeatureMapper {

    public MNPMapper(FlowFeatureMapperArgumentCollection fmArgs) {
        super(fmArgs);
    }

    @Override
    protected FlowFeatureMapper.MappedFeature detectFeature(GATKRead read, ReferenceContext referenceContext, int readOfs, int refOfs) {
        final byte[] bases = read.getBasesNoCopy();
        final byte[] ref = referenceContext.getBases();

        if ( ref[refOfs] != 'N' && (bases[readOfs] != ref[refOfs]) ) {

            // SNV located, next scan for a second SNV within the allowed distance
            int featureSize = 2;
            int readOfs2 = readOfs + 1;
            int refOfs2 = refOfs + 1;
            boolean found = false;
            for ( ; featureSize <= fmArgs.maxMnpDistance + 2 ; featureSize++, readOfs2++, refOfs2++ ) {

                // still within read/ref limits
                if ( readOfs2  >= bases.length || refOfs2 >= ref.length ) {
                    return null;
                }

                // if found difference, break out of loop
                if ( ref[refOfs2] == 'N' ) {
                    return null;
                }
                if ( bases[readOfs2] != ref[refOfs2] ) {
                    found = true;
                    break;
                }
            }
            if ( !found ) {
                return null;
            }

            // check if surrounded
            boolean surrounded = isSurrounded(read, referenceContext, readOfs, refOfs, featureSize);
            if ( ignoreBecauseNotSurrounded(surrounded) ) {
                return null;
            }

            // add this feature
            byte[] readBases = Arrays.copyOfRange(bases, readOfs, readOfs + featureSize);
            byte[] refBases = Arrays.copyOfRange(ref, refOfs, refOfs + featureSize);

            FlowFeatureMapper.MappedFeature feature = FlowFeatureMapper.MappedFeature.makeFeature(
                    FlowFeatureMapperArgumentCollection.MappingFeatureEnum.SNV,
                    read, readBases, refBases, readOfs, referenceContext.getStart() + refOfs, readOfs - refOfs);

            return enrichFeature(feature, surrounded);
        } else {
            return null;
        }
    }
}
