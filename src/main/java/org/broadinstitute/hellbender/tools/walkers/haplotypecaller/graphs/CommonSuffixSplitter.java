package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.utils.Utils;

//import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

/**
 * Split a collection of middle nodes in a graph into their shared prefix and suffix values
 *
 * This code performs the following transformation.  Suppose I have a set of vertices V, such
 * that each vertex is composed of sequence such that
 *
 * Vi = prefix + seq_i + suffix
 *
 * where prefix and suffix are shared sequences across all vertices V.  This replaces each
 * Vi with three nodes prefix, seq_i, and suffix connected in a simple chain.
 *
 * This operation can be performed in a very general case, without too much worry about the incoming
 * and outgoing edge structure of each Vi.  The partner algorithm SharedSequenceMerger can
 * put these pieces back together in a smart way that maximizes the sharing of nodes
 * while respecting complex connectivity.
 *
 */
public final class CommonSuffixSplitter {

    private CommonSuffixSplitter() {}

    /**
     * Simple single-function interface to split and then update a graph
     *
     * @param graph the graph containing the vertices in toMerge
     * @param v The bottom node whose incoming vertices we'd like to split
     * @return true if some useful splitting was done, false otherwise
     */
    public static boolean split(final SeqGraph graph, final SeqVertex v) {
        Utils.nonNull(graph, "graph cannot be null");
        Utils.nonNull(v, "v cannot be null");
        Utils.validateArg(graph.vertexSet().contains(v), () -> "graph doesn't contain vertex v " + v);
        final Collection<SeqVertex> toSplit= graph.incomingVerticesOf(v);

        final SeqVertex suffixVTemplate = commonSuffix(graph, v, toSplit, true);
        if (suffixVTemplate == null){
            return false;
        }
        final List<BaseEdge> edgesToRemove = new LinkedList<>();

//                graph.printGraph(new File("split.pre_" + v.getSequenceString() + "." + counter + ".dot"), 0);
                for ( final SeqVertex mid : toSplit ) {
                    // create my own copy of the suffix
                    final SeqVertex suffixV = new SeqVertex(suffixVTemplate.getSequence());
                    graph.addVertex(suffixV);
                    final SeqVertex prefixV = mid.withoutSuffix(suffixV.getSequence());
                    final BaseEdge out = graph.outgoingEdgeOf(mid);

                    final SeqVertex incomingTarget;
                    if ( prefixV == null ) {
                        // this node is entirely explained by suffix
                        incomingTarget = suffixV;
                    } else {
                        incomingTarget = prefixV;
                        graph.addVertex(prefixV);
                        graph.addEdge(prefixV, suffixV, new BaseEdge(out.isRef(), 1,1));
                        edgesToRemove.add(out);
                    }

                    graph.addEdge(suffixV, graph.getEdgeTarget(out), out.copy());

                    for ( final BaseEdge in : graph.incomingEdgesOf(mid) ) {
                        graph.addEdge(graph.getEdgeSource(in), incomingTarget, in.copy());
                        edgesToRemove.add(in);
                    }
                }

                graph.removeAllVertices(toSplit);
                graph.removeAllEdges(edgesToRemove);
//                graph.printGraph(new File("split.post_" + v.getSequenceString() + "." + counter++ + ".dot"), 0);

                return true;
    }

    /**
     * Simple single-function interface to split and then update a graph
     *
     * @param graph the graph containing the vertices in toSplit
     * @param v The top node whose outgoing vertices we'd like to split
     * @return true if some useful splitting was done, false otherwise
     */
    public static boolean splitFromTop(final SeqGraph graph, final SeqVertex v) {
        Utils.nonNull(graph, "graph cannot be null");
        Utils.nonNull(v, "v cannot be null");
        Utils.validateArg(graph.vertexSet().contains(v), () -> "graph doesn't contain vertex v " + v);
        final Collection<SeqVertex> toSplit= graph.outgoingVerticesOf(v);

        final SeqVertex suffixVTemplate = commonSuffix(graph, v, toSplit,false);
        if (suffixVTemplate == null){
            return false;
        }
        final List<BaseEdge> edgesToRemove = new LinkedList<>();
        final List<SeqVertex> nodesToMerge = new LinkedList<>();

//        graph.printGraph(new File("topsplit.pre_" + v.getSequenceString() + "." + counter + ".dot"), 0);
        for ( final SeqVertex mid : toSplit ) {
            // create my own copy of the suffix
            final SeqVertex suffixV = new SeqVertex(suffixVTemplate.getSequence());
            nodesToMerge.add(suffixV);
            graph.addVertex(suffixV);
            final SeqVertex prefixV = mid.withoutSuffix(suffixV.getSequence());
            final BaseEdge in = graph.incomingEdgeOf(mid);
            final SeqVertex topTarget;
            if ( prefixV == null ) {
                // this node is entirely explained by suffix
                topTarget=suffixV;
            } else {
                topTarget=prefixV;
                graph.addVertex(prefixV);
                graph.addEdge(prefixV, suffixV, in.copy());
                edgesToRemove.add(in);
            }
            graph.addEdge(v, topTarget, in.copy());

            graph.addEdge(suffixV, graph.getEdgeTarget(in), in.copy());

            for ( final BaseEdge out : graph.outgoingEdgesOf(mid) ) {
                graph.addEdge(suffixV,graph.getEdgeTarget(out), out.copy());
                edgesToRemove.add(out);
            }
        }

        graph.removeAllVertices(toSplit);
        graph.removeAllEdges(edgesToRemove);
//        graph.printGraph(new File("topsplit.mid_" + v.getSequenceString() + "." + counter + ".dot"), 0);

        merge(graph,nodesToMerge);

//        graph.printGraph(new File("topsplit.post_" + v.getSequenceString() + "." + counter++ + ".dot"), 0);

        return true;
    }


    private static void merge(final SeqGraph graph, final List<SeqVertex> nodesToMerge) {
        final Collection<BaseEdge> edgesToRemove = new LinkedList<>();

        final SeqVertex newVertex = new SeqVertex(nodesToMerge.get(0).getSequence());
        graph.addVertices(newVertex);

        for (final SeqVertex node : nodesToMerge) {
            for (final BaseEdge edge : new ArrayList<>(graph.incomingEdgesOf(node))) {
                GraphUtils.addOrModifyEdge(graph,graph.getEdgeSource(edge),newVertex,edge);
                edgesToRemove.add(edge);
            }
        }

        for (final SeqVertex node : nodesToMerge) {
            for (final BaseEdge edge :  new ArrayList<>(graph.outgoingEdgesOf(node))) {
                edgesToRemove.add(edge);
                GraphUtils.addOrModifyEdge(graph, newVertex, graph.getEdgeTarget(edge), edge);
            }
        }
        graph.removeAllEdges(edgesToRemove);
        graph.removeAllVertices(nodesToMerge);
    }

    private static SeqVertex commonSuffix(final SeqGraph graph, final SeqVertex v, final Collection<SeqVertex> toSplit, final boolean bottom) {
        if (toSplit.size() < 2) {
            // Can only split at least 2 vertices
            return null;
        } else if (!(bottom ?
                safeToSplitBot(graph, v, toSplit) :
                safeToSplitTop(graph, v, toSplit))) {
            return null;
        }

        final SeqVertex suffixVTemplate = commonSuffix(toSplit);
        if ( suffixVTemplate.isEmpty() ) {
            return null;
        } else if ( wouldEliminateRefSource(graph, suffixVTemplate, toSplit) ) {
            return null;
        } else if ( allVerticesAreTheCommonSuffix(suffixVTemplate, toSplit) ) {
            return null;
        } else {
            return suffixVTemplate;
        }
    }

    /**
     * Would factoring out this suffix result in eliminating the reference source vertex?
     * @param graph the graph
     * @param commonSuffix the common suffix of all toSplits
     * @param toSplits the list of vertices we're are trying to split
     * @return true if toSplit contains the reference source and this ref source has all and only the bases of commonSuffix
     */
    private static boolean wouldEliminateRefSource(final SeqGraph graph, final SeqVertex commonSuffix, final Collection<SeqVertex> toSplits) {
        for ( final SeqVertex toSplit : toSplits ) {
            if ( graph.isRefSource(toSplit) ) {
                return toSplit.length() == commonSuffix.length();
            }
        }
        return false;
    }

//    private static int counter = 0;

    /**
     * Would all vertices that we'd split just result in the common suffix?
     *
     * That is, suppose we have prefix nodes ABC and ABC.  After splitting all of the vertices would
     * just be ABC again, and we'd enter into an infinite loop.
     *
     * @param commonSuffix the common suffix of all vertices in toSplits
     * @param toSplits the collection of vertices we want to split
     * @return true if all of the vertices are equal to the common suffix
     */
    private static boolean allVerticesAreTheCommonSuffix(final SeqVertex commonSuffix, final Collection<SeqVertex> toSplits) {
        for ( final SeqVertex toSplit : toSplits ) {
            if ( toSplit.length() != commonSuffix.length() ) {
                return false;
            }
        }

        return true;
    }

    /**
     * Can we safely split up the vertices in toMerge?
     *
     * @param graph a graph
     * @param bot a vertex whose incoming vertices we want to split
     * @param toMerge the set of vertices we'd be splitting up
     * @return true if we can safely split up toMerge
     */
    private static boolean safeToSplitBot(final SeqGraph graph, final SeqVertex bot, final Collection<SeqVertex> toMerge) {
        final Collection<SeqVertex> outgoingOfBot = new HashSet<>(graph.outgoingVerticesOf(bot));
        for ( final SeqVertex m : toMerge ) {

            final Set<BaseEdge> outs = graph.outgoingEdgesOf(m);
            if ( m == bot || outs.size() != 1 || ! graph.outgoingVerticesOf(m).contains(bot) )
                // m == bot => don't allow self cycles in the graph
            {
                return false;
            }
            if ( outgoingOfBot.contains(m) )
                // forbid cycles from bottom -> mid
            {
                return false;
            }
        }

        return true;
    }
    /**
     * Can we safely split up the vertices in toMerge?
     *
     * @param graph a graph
     * @param top a vertex whose incoming vertices we want to split
     * @param toMerge the set of vertices we'd be splitting up
     * @return true if we can safely split up toMerge
     */
    private static boolean safeToSplitTop(final SeqGraph graph, final SeqVertex top, final Collection<SeqVertex> toMerge) {
        final Collection<SeqVertex> incommingOfTop = new HashSet<>(graph.incomingVerticesOf(top));
        for (final SeqVertex m : toMerge) {

            final Set<BaseEdge> ins = graph.incomingEdgesOf(m);
            if (m == top || ins.size() != 1 || !graph.incomingVerticesOf(m).contains(top)) {
                // m == top => don't allow self cycles in the graph
                return false;
            }
            if (incommingOfTop.contains(m)) {
                // forbid cycles from bottom -> mid
                return false;
            }
        }

        return true;
    }


    /**
     * Return the longest suffix of bases shared among all provided vertices
     *
     * For example, if the vertices have sequences AC, CC, and ATC, this would return
     * a single C.  However, for ACC and TCC this would return CC.  And for AC and TG this
     * would return null;
     *
     * @param middleVertices a non-empty set of vertices
     * @return a single vertex that contains the common suffix of all middle vertices
     */
    private static SeqVertex commonSuffix(final Collection<SeqVertex> middleVertices) {
        final List<byte[]> kmers = GraphUtils.getKmers(middleVertices);
        final int min = GraphUtils.minKmerLength(kmers);
        final int suffixLen = GraphUtils.commonMaximumSuffixLength(kmers, min);
        final byte[] kmer = kmers.get(0);
        final byte[] suffix = Arrays.copyOfRange(kmer, kmer.length - suffixLen, kmer.length);
        return new SeqVertex(suffix);
    }
}