package org.ultimagenomics.variant_calling;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceMode;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;

import java.util.Collection;
import java.util.List;


/**
 * Calculate likelyhood matrix for each Allele in VCF against a set of Reads limited by a set of Haplotypes
 *
 * <h3>How HaplotypeBasedVariantRecaller works</h3>
 * For every variant in the VCF:
 * <ul>
 *     <li>Fetch all haplotypes that span it </li>
 *     <li>Trim the reads to the haplotypes </li>
 *     <li>Calculate ReadLikelihoodMatrix as usual in HaplotypeCallerEngine. Do not apply filterPoorlyModeledReads in process </li>
 *     <li>assignGenotypeLikelihoods of HaplotypeCallerGenotypingEngine, simillar to HaplotypeCaller</li>
 *     <li>Print this output </li>
 * </ul>
 * <br />
 *
 * <h3>Input</h3>
 * <ul>
 *     <li>Input VCG file with Alelles to recall</li>
 *     <li>Input BAM file with Reads against which Alelles are recalled</li>
 *     <li>Input BAM file with Haplotypes to limit reads by</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <p>
 *     Output matrix file (multipart csv)
 * </p>
 *
 */
@CommandLineProgramProperties(
        summary = "Recalling variants from haplotypes",
        oneLineSummary = "Calculate likelyhood matrix for each Allele in VCF against a set of Reads limited by a set of Haplotypes",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public final class HaplotypeBasedVariantRecaller extends GATKTool {

    private static final Logger logger = LogManager.getLogger(HaplotypeBasedVariantRecaller.class);

    @ArgumentCollection
    private HaplotypeBasedVariantRecallerArgumentCollection vrArgs = new HaplotypeBasedVariantRecallerArgumentCollection();

    @Override
    public void traverse() {
        logger.info("traverse() called");
    }
}
