package org.ultimagenomics.variant_recalling;

import org.broadinstitute.barclay.argparser.Argument;

import java.io.File;
import java.io.Serializable;

/**
 * Set of arguments for the {@link HaplotypeBasedVariantRecaller}
 */
public class HaplotypeBasedVariantRecallerArgumentCollection implements Serializable{
    private static final long serialVersionUID = 1L;

    /**
     *  This argument specifies a VCF file with Alleles to be recalled
     **/
    @Argument(fullName = "alleles-file-vcf", doc = "VCF file containing alleles", optional = false)
    public String ALLELE_VCF_FILE = null;

    /**
     *  This argument specifies a BAM file with Haplotypes to limit reads on
     **/
    @Argument(fullName = "haplotypes-file-bam", doc = "BAM file containing haplotypes", optional = false)
    public String HAPLOTYPES_BAM_FILE = null;

    /**
     *  This argument specifies a CSV to be filled with likelihood matrix data
     **/
    @Argument(fullName = "matrix-file-csv", doc = "CSV file to be filled with likelihood matrix data", optional = false)
    public File MATRIX_CSV_FILE = null;
}
