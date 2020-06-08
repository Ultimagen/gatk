package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Collection;
import java.util.Comparator;

/**
 * Simple edge class for connecting nodes in the graph.
 *
 * Works equally well for all graph types (kmer or sequence)
 */
public class BaseEdge {
    protected int forwardMultiplicity;
    protected int reverseMultiplicity;
    private boolean isRef;

    /**
     * Create a new BaseEdge with weight multiplicity and, if isRef == true, indicates a path through the reference
     *
     * @param isRef indicates whether this edge is a path through the reference
     * @param forwardMultiplicity the number of observations of this edge
     * @param reverseMultiplicity the number of observations of this edge
     */
    public BaseEdge(final boolean isRef, final int forwardMultiplicity, final int reverseMultiplicity) {
        Utils.validateArg( forwardMultiplicity >= 0, () -> "multiplicity must be >= 0 but got " + forwardMultiplicity);
        Utils.validateArg( reverseMultiplicity >= 0, () -> "multiplicity must be >= 0 but got " + reverseMultiplicity);
        this.forwardMultiplicity = forwardMultiplicity;
        this.reverseMultiplicity = reverseMultiplicity;
        this.isRef = isRef;
    }

    /**
     * Create a new copy of this BaseEdge
     */
    public BaseEdge copy() {
        return new BaseEdge(isRef(), forwardMultiplicity, reverseMultiplicity);
    }

    /**
     * Get the number of observations of paths connecting two vertices
     * @return a positive integer >= 0
     */
    public final int getMultiplicity() {
        return forwardMultiplicity + reverseMultiplicity;
    }

    /**
     * Get the DOT format label for this edge, to be displayed when printing this edge to a DOT file
     * @return a non-null string
     */
    public String getDotLabel() {
        return String.format("%d//%d",forwardMultiplicity, reverseMultiplicity);
    }

    /**
     * Increase the multiplicity of this edge by incr
     * @param incr the change in the multiplicity, must be >= 0
     * @param reverseStrand whether the change in the multiplicity shoudl be attributed to the forward or reverse strand
     */
    public void incMultiplicity(final int incr, final boolean reverseStrand) {
        Utils.validateArg(incr >= 0, () -> "incr must be >= 0 but got " + incr);
        this.forwardMultiplicity += reverseStrand ? 0 : incr;
        this.reverseMultiplicity += reverseStrand ? incr : 0;

    }

    public void incMultiplicity(final BaseEdge other) {
        this.forwardMultiplicity += other.forwardMultiplicity;
        this.reverseMultiplicity += other.reverseMultiplicity;
    }

    /**
     * A special accessor that returns the multiplicity that should be used by pruning algorithm
     *
     * @return the multiplicity value that should be used for pruning
     */
    public int getPruningMultiplicity() {
            return getMultiplicity();
    }

    /**
     * Set the multiplicity of this edge to value
     * @param forwardMultiplicity an integer >= 0
     * @param reverseMultiplicity an integer >= 0
     */
    public final void setMultiplicity( final int forwardMultiplicity, final int reverseMultiplicity) {
        ParamUtils.isPositiveOrZero(forwardMultiplicity, "forwardMultiplicity must be >= 0, but got " + forwardMultiplicity);
        ParamUtils.isPositiveOrZero(reverseMultiplicity, "reverseMultiplicity must be >= 0, but got " + reverseMultiplicity);
        this.forwardMultiplicity = forwardMultiplicity;
        this.reverseMultiplicity = reverseMultiplicity;
    }

    public final void setMultiplicity(final BaseEdge other){
        setMultiplicity(other.forwardMultiplicity,other.reverseMultiplicity);
    }

    /**
     * Does this edge indicate a path through the reference graph?
     * @return true if so
     */
    public final boolean isRef() {
        return isRef;
    }

    /**
     * Indicate that this edge follows the reference sequence, or not
     * @param isRef true if this is a reference edge
     */
    public final void setIsRef( final boolean isRef ) {
        this.isRef = isRef;
    }

    /**
     * Sorts a collection of BaseEdges in decreasing order of weight, so that the most
     * heavily weighted is at the start of the list
     */
    //TODO: this doesn't seem to be used anywhere....
    public static final Comparator<BaseEdge> EDGE_MULTIPLICITY_ORDER = Comparator.comparingInt(BaseEdge::getMultiplicity).reversed();

    /**
     * Add edge to this edge, updating isRef and multiplicity as appropriate
     *
     * isRef is simply the or of this and edge
     * multiplicity is the sum
     *
     * @param edge the edge to add
     * @return this
     */
    public final BaseEdge add(final BaseEdge edge) {
        Utils.nonNull(edge, "edge cannot be null");
        this.forwardMultiplicity += edge.forwardMultiplicity;
        this.reverseMultiplicity += edge.reverseMultiplicity;
        isRef = isRef || edge.isRef();
        return this;
    }

    /**
     * Create a new BaseEdge with the multiplicity equal to the sum of the multiplicities. The resulting edge is a reference edge if any of the argument edges are reference.
     *
     * @param edges a collection of edges to or their isRef values
     * @return a newly allocated BaseEdge
     */
    public static BaseEdge makeOREdge(final Collection<BaseEdge> edges) {
        Utils.nonNull(edges);
        final boolean anyRef = edges.stream().anyMatch(BaseEdge::isRef);
        final int forwardSum = edges.stream().mapToInt(e -> e.forwardMultiplicity).sum();
        final int reverseSum = edges.stream().mapToInt(e -> e.reverseMultiplicity).sum();

        return new BaseEdge(anyRef, forwardSum, reverseSum);
    }

    @Override
    public final String toString() {
        return String.format("BaseEdge{multiplicity=%d//%d, isRef=%b}", forwardMultiplicity, reverseMultiplicity, isRef);
    }
}
