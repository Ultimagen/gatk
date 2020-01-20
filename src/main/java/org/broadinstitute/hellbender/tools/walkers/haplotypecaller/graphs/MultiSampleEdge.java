package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.PriorityQueue;

/**
 * Edge class for connecting nodes in the graph that tracks some per-sample information.
 * <p>
 * This class extends BaseEdge with the additional functionality of tracking the maximum
 * multiplicity seen within any single sample.  The workflow for using this class is:
 * </p>
 * <pre>
 * {@code
 *      MultiSampleEdge e = new MultiSampleEdge(ref, 1)
 *      e.incMultiplicity(1)              // total is 2, per sample is 2, max per sample is 1
 *      e.getPruningMultiplicity()        // = 1
 *      e.flushSingleSampleMultiplicity() // total is 2, per sample is 0, max per sample is 2
 *      e.getPruningMultiplicity()        // = 2
 *      e.incMultiplicity(3)              // total is 5, per sample is 3, max per sample is 2
 *      e.getPruningMultiplicity()        // = 2
 *      e.flushSingleSampleMultiplicity() // total is 5, per sample is 0, max per sample is 3
 *      e.getPruningMultiplicity()        // = 3
 * }
 * </pre>
 */
public final class MultiSampleEdge extends BaseEdge {
    private BaseEdge currentSingleSampleEdge;
    private final int singleSampleCapacity;
    private final PriorityQueue<BaseEdge> singleSampleMultiplicities;

    private final List<Integer> referencePathIndexes = new ArrayList<>(2);

    /**
     * Create a new MultiSampleEdge with weight multiplicity and, if isRef == true, indicates a path through the reference
     *
     * @param isRef                indicates whether this edge is a path through the reference
     * @param forwardMult          the number of forward observations of this edge
     * @param reverseMult          the number of reverse observations of this edge
     * @param singleSampleCapacity the max number of samples to track edge multiplicities
     */
    public MultiSampleEdge(final boolean isRef, final int forwardMult, final int reverseMult, final int singleSampleCapacity) {
        super(isRef, forwardMult, reverseMult);

        Utils.validateArg( singleSampleCapacity > 0, () -> "singleSampleCapacity must be > 0 but found: " + singleSampleCapacity);
        singleSampleMultiplicities = new PriorityQueue<>(singleSampleCapacity, new BaseEdgeComparator());
        singleSampleMultiplicities.add(new BaseEdge(false, forwardMult,reverseMult));

        currentSingleSampleEdge = new BaseEdge(isRef, forwardMult,reverseMult);

        this.singleSampleCapacity = singleSampleCapacity;
    }

    private static class BaseEdgeComparator implements Comparator<BaseEdge> {
        @Override
        public int compare(final BaseEdge o1, final BaseEdge o2) {
            return o1.getPruningMultiplicity() - o2.getPruningMultiplicity();
        }
    }

    @Override
    public MultiSampleEdge copy() {
        return new MultiSampleEdge(isRef(), forwardMultiplicity, reverseMultiplicity,singleSampleCapacity); // TODO -- should I copy values for other features?
    }

    /**
     * update the single sample multiplicities by adding the current single sample multiplicity to the priority queue, and
     * reset the current single sample multiplicity to 0.
     */
    public void flushSingleSampleMultiplicity() {
        singleSampleMultiplicities.add(currentSingleSampleEdge);
        if (singleSampleMultiplicities.size() == singleSampleCapacity + 1) {
            singleSampleMultiplicities.poll(); // remove the lowest multiplicity from the list
        } else if (singleSampleMultiplicities.size() > singleSampleCapacity + 1) {
            throw new IllegalStateException("Somehow the per sample multiplicity list has grown too big: " + singleSampleMultiplicities);
        }
        currentSingleSampleEdge = new BaseEdge(isRef(), 0, 0);
    }

    @Override
    public void incMultiplicity(final int incr, final boolean fromReverseStrand) {
        super.incMultiplicity(incr, fromReverseStrand);
        currentSingleSampleEdge.incMultiplicity(incr,fromReverseStrand);
    }

    @Override
    public int getPruningMultiplicity() {
        return currentSingleSampleEdge.getPruningMultiplicity() +
                Optional.ofNullable(singleSampleMultiplicities.peek())
                        .orElse(new BaseEdge(isRef(), 0, 0))
                        .getPruningMultiplicity();
    }

    @Override
    public String getDotLabel() {
        return super.getDotLabel() + "->" + getPruningMultiplicity();
    }

    @VisibleForTesting
    int getCurrentSingleSampleMultiplicity() {
        return currentSingleSampleEdge.getMultiplicity();
    }

    public void addReferenceIndex(int i) {
        referencePathIndexes.add(i);
    }

    public List<Integer> getReferencePathIndexes() {
        return referencePathIndexes;
    }
}
