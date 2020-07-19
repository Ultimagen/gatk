package org.ultimagenomics.variant_calling;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.genotyper.StandardCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerReadThreadingAssemblerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadThreadingAssemblerArgumentCollection;
import org.ultimagenomics.flow_based_read.utils.FlowBasedAlignmentArgumentCollection;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Set of arguments for the {@link HaplotypeBasedVariantRecaller}
 */
public class HaplotypeBasedVariantRecallerArgumentCollection implements Serializable{
    private static final long serialVersionUID = 1L;

    /**
     *  This argument specifies a VCF file with Alleles to be recalled
     **/
    @Argument(fullName = "alleles-file-vcf", doc = "VCF file containing alleles", optional = false)
    public File ALLELE_VCF_FILE = null;

    /**
     *  This argument specifies a BAM file with Reads to recall from
     **/
    @Argument(fullName = "reads-file-bam", doc = "BAM file containing reads", optional = false)
    public File READS_BAM_FILE = null;

    /**
     *  This argument specifies a BAM file with Haplotypes to limit reads on
     **/
    @Argument(fullName = "haplotypes-file-bam", doc = "BAM file containing haplotypes", optional = false)
    public File HAPLOTYPES_BAM_FILE = null;

    /**
     *  This argument specifies a CSV to be filled with likelihood matrix data
     **/
    @Argument(fullName = "matrix-file-csv", doc = "CSV file to be filled with likelihood matrix data", optional = false)
    public File MATRIX_CSV_FILE = null;

    /**
     *  This argument specifies FASTA reference
     **/
    @Argument(fullName = "reference-fasta", doc = "FASTA reference", optional = false)
    public File REFERENCE_FASTA = null;

    /**
     *  This argument specifies region to work on
     **/
    @Argument(fullName = "region-loc", doc = "Region to work on", optional = false)
    public String REGION_LOC = null;
}
