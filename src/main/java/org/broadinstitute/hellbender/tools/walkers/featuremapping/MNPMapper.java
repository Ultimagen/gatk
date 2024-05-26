package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.function.Consumer;

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
    public void forEachOnRead(GATKRead read, ReferenceContext referenceContext, Consumer<? super FlowFeatureMapper.MappedFeature> action) {

    }

    public FilterStatus noFeatureButFilterAt(GATKRead read, ReferenceContext referenceContext, int start) {

        return FilterStatus.None;
    }

}
