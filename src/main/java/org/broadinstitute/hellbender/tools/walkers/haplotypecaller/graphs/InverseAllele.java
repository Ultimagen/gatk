package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import htsjdk.variant.variantcontext.Allele;

public class InverseAllele extends Allele {
    final static public long serialVersionUID = 1L;

    private final Allele internalAllele;

    private InverseAllele(final Allele allele) {
        super(allele, false);
        this.internalAllele = allele;
    }

    public static Allele of(final Allele allele){
        if (allele instanceof InverseAllele) {
            return ((InverseAllele)allele).internalAllele;
        } else {
            return new InverseAllele(allele);
        }
    }

    @Override
    public boolean isReference() {
        return internalAllele.isReference();
    }

    @Override
    public boolean isSymbolic() {
        return true;
    }

    @Override
    public String getDisplayString() {
        return "~" + super.getDisplayString();
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof InverseAllele) )
            return false;

        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        if (!super.equals(o)) {
            return false;
        }

        final InverseAllele that = (InverseAllele) o;

        return internalAllele.equals(that.internalAllele);
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + internalAllele.hashCode();
        return result;
    }
}
