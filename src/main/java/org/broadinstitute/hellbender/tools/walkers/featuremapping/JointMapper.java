package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.function.Consumer;

/**
 * An implementation of a feature mapper that aggregates several other feature mappers
 *
 */

public class JointMapper implements FeatureMapper {

    final List<FeatureMapper> fmList;

    public JointMapper(List<FeatureMapper> fmList) {
        this.fmList = fmList;
    }

    @Override
    public void forEachOnRead(GATKRead read, ReferenceContext referenceContext, Consumer<? super FlowFeatureMapper.MappedFeature> action) {
        for ( FeatureMapper fm : fmList ) {
            fm.forEachOnRead(read, referenceContext, action);
        }
    }

    public FilterStatus noFeatureButFilterAt(GATKRead read, ReferenceContext referenceContext, int start) {

        for ( FeatureMapper fm : fmList ) {
            FilterStatus status = fm.noFeatureButFilterAt(read, referenceContext, start);
            if ( status != FilterStatus.None ) {
                return status;
            }
        }

        return FilterStatus.None;
    }

}
