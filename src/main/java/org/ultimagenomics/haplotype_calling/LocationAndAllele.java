package org.ultimagenomics.haplotype_calling;
import htsjdk.variant.variantcontext.Allele;

public class LocationAndAllele extends Allele {
    final static public long serialVersionUID = 1L;
    private final int loc;

    public LocationAndAllele(final int loc, final Allele allele) {
        super(allele, false);
        this.loc = loc;

    }

    public int getLoc() {
        return loc;
    }

    public Allele getAllele() {
        return this;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final org.ultimagenomics.haplotype_calling.LocationAndAllele that = (org.ultimagenomics.haplotype_calling.LocationAndAllele) o;

        if (loc != that.loc) return false;
        return this != null ? super.equals(that) : that == null;
    }

    @Override
    public int hashCode() {
        return 31 * loc + (this != null ? super.hashCode() : 0);
    }

    public String toString() {return String.format("(%d) %s", loc, getBaseString());}
}

