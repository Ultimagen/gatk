package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import htsjdk.variant.variantcontext.Allele;

public class BaseVertexBackedAllele<T extends BaseVertex> extends Allele {
    final static public long serialVersionUID = 1L;

    final private T baseVertex;

    static private byte[] NULL_ALLELE_BYTES = ".".getBytes();

    public BaseVertexBackedAllele(T vertex) {
        super(vertex.sequence.length == 0 ? NULL_ALLELE_BYTES : vertex.sequence, false);
        this.baseVertex = vertex;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }

        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        @SuppressWarnings("unchecked")
        final BaseVertexBackedAllele<T> that = (BaseVertexBackedAllele<T>) o;

        return baseVertex.equals(that.baseVertex);
    }

    @Override
    public int hashCode() {
        return baseVertex.hashCode();
    }
}
