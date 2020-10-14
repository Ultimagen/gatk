package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import htsjdk.variant.variantcontext.Allele;

public class JoinedContigs extends Allele {
    final static public long serialVersionUID = 1L;


    private final Allele allele1;
    private final Allele allele2;

    public JoinedContigs(final Allele allele1, final Allele allele2) {
        super("<"+allele1.getBaseString() + "->" + allele2.getBaseString()+">", false);
        this.allele1 = allele1;
        this.allele2 = allele2;
    }

    @Override
    public boolean isReference() {
        return allele1.isReference() && allele2.isReference();
    }

    @Override
    public byte [] getBases(){
        return super.getBases();
    }
    @Override
    public boolean isSymbolic() {
        return true;
    }

    public Allele getAllele1() {
        return allele1;
    }

    public Allele getAllele2() {
        return allele2;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        if (!super.equals(o)) {
            return false;
        }

        final JoinedContigs that = (JoinedContigs) o;

        if (!that.allele1.equals(this.allele1)) {
            return false;
        }
        return that.allele2.equals(this.allele2);
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + allele1.hashCode();
        result = 31 * result + allele2.hashCode();
        return result;
    }
}
