package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * An implementation of a feature mapper that finds SNPs (SVN)
 *
 * This class only finds SNP that are surrounded by a specific number of bases identical to the reference.
 */

public class SNVMapper extends BaseFeatureMapper implements FeatureMapper {

    public SNVMapper(FlowFeatureMapperArgumentCollection fmArgs, SAMFileHeader hdr) {
        super(fmArgs, hdr);
    }

    @Override
    protected FlowFeatureMapper.MappedFeature detectFeature(GATKRead read, ReferenceContext referenceContext, final byte bases[], final byte ref[], int readOfs, int refOfs) {

        if ( ref[refOfs] != 'N' && (fmArgs.reportAllAlts || (bases[readOfs] != ref[refOfs])) ) {

            // check if surrounded
            boolean surrounded = isSurrounded(read, referenceContext, readOfs, refOfs, 1, 1);
            if ( ignoreBecauseNotSurrounded(surrounded) ) {
                return null;
            }

            // add this feature
            FlowFeatureMapper.MappedFeature feature = FlowFeatureMapper.MappedFeature.makeFeature(
                    FlowFeatureMapperArgumentCollection.MappingFeatureEnum.SNV,
                        read, readOfs, ref[refOfs], referenceContext.getStart() + refOfs, readOfs - refOfs);

            return enrichFeature(feature, surrounded);
        } else {
            return null;
        }
    }

}
