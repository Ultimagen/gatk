package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.Utils;

public final class ReadDirectionAwareEdge extends BaseEdge {
    private int forwardMultiplicity;
    private int reverseMultiplicity;

    /**
     * Create a new MultiSampleEdge with weight multiplicity and, if isRef == true, indicates a path through the reference
     *
     * @param isRef indicates whether this edge is a path through the reference
     */
    public ReadDirectionAwareEdge(final boolean isRef, final int forwardMultiplicity, final int reverseMultiplicity) {
        super(isRef, forwardMultiplicity, reverseMultiplicity);
        Utils.validateArg(forwardMultiplicity > 0, () -> "forwardMultiplicity must be > 0 but found: " + forwardMultiplicity);
        Utils.validateArg(reverseMultiplicity > 0, () -> "reverseMultiplicity must be > 0 but found: " + reverseMultiplicity);
        this.forwardMultiplicity = forwardMultiplicity;
        this.reverseMultiplicity = reverseMultiplicity;
    }

    @Override
    public ReadDirectionAwareEdge copy() {
        return new ReadDirectionAwareEdge(isRef(), forwardMultiplicity, reverseMultiplicity);
    }
}
