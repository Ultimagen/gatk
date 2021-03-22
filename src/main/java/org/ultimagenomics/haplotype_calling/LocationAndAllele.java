package org.ultimagenomics.haplotype_calling;
import htsjdk.variant.variantcontext.Allele;

/**
 * This class is similar to {@link org.broadinstitute.hellbender.tools.walkers.haplotypecaller.LocationAndAlleles} but
 * allows to keep only an allele/ref pair rather than a list of alleles. The comparison is done on allele by allele basis and
 * not in the way it is done on LocationAndAlleles
 */

public class LocationAndAllele extends Allele {
    final static public long serialVersionUID = 1L;
    private final int loc;
    private final String contig;
    private final Allele refAllele;
    public LocationAndAllele(final String contig, final int loc, final Allele allele, final Allele refAllele) {
        super(allele, false);
        this.loc = loc;
        this.contig = contig;
        this.refAllele = refAllele;
    }

    public int getLoc() {
        return loc;
    }

    public String getContig() { return contig; }

    public Allele getAllele() {
        return this;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final org.ultimagenomics.haplotype_calling.LocationAndAllele that = (org.ultimagenomics.haplotype_calling.LocationAndAllele) o;

        if (loc != that.loc) return false;
        return super.equals(that) && this.refAllele.equals(that.getRefAllele());
    }

    @Override
    public int hashCode() {
        return 31 * loc + (this != null ? super.hashCode() : 0);
    }

    public String toString() {return String.format("(%d) %s/%s", loc, getBaseString(), getRefAllele().getBaseString());}
    public Allele getRefAllele() { return refAllele;}
}

