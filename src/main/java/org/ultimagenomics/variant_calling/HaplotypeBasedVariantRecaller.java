package org.ultimagenomics.variant_calling;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.lang3.tuple.Pair;
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
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.ultimagenomics.flow_based_read.alignment.FlowBasedAlignmentEngine;
import org.ultimagenomics.flow_based_read.read.FlowBasedRead;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;


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

    @ArgumentCollection
    private final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();


    private ReadLikelihoodCalculationEngine likelihoodCalculationEngine = null;

    @Override
    public void traverse() {
        logger.info("traverse() called");
        logger.info("Alleles VCF: " + vrArgs.ALLELE_VCF_FILE);
        logger.info("Reads BAM: " + vrArgs.READS_BAM_FILE);
        logger.info("Haplotypes BAM: " + vrArgs.HAPLOTYPES_BAM_FILE);
        logger.info("Matrix CSV: " + vrArgs.MATRIX_CSV_FILE);

        // inits
        likelihoodCalculationEngine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(hcArgs.likelihoodArgs, hcArgs.fbargs);
        String[] sampleNames = {"sm1"};
        SampleList samplesList = new IndexedSampleList(Arrays.asList(sampleNames));
        HaplotypeCallerGenotypingEngine genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samplesList, !hcArgs.doNotRunPhysicalPhasing);
        ReferenceDataSource             reference = ReferenceDataSource.of(vrArgs.REFERENCE_FASTA.toPath());
        VariantCallerResultWriter       resultWriter = new VariantCallerResultWriter(vrArgs.MATRIX_CSV_FILE);

        // walk regions defined by haploype groups
        final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<VariantContext>(
                vrArgs.ALLELE_VCF_FILE.getAbsolutePath(), null, 0, VariantContext.class);
        final HaplotypeRegionWalker regionWalker = new HaplotypeRegionWalker(vrArgs);
        final TrimmedReadsReader readsReader = new TrimmedReadsReader(vrArgs);
        regionWalker.forEach(haplotypes -> {

            // get reads overlapping haplotype, get variants
            SimpleInterval  loc = new SimpleInterval(haplotypes.get(0).getGenomeLocation());
            Collection<FlowBasedRead> reads = readsReader.getReads(loc);
            List<VariantContext> variants = dataSource.queryAndPrefetch(loc);
            logger.info(String.format("%s: %d haplotypes, %d reads, %d variants",
                    loc.toString(), haplotypes.size(), reads.size(), variants.size()));

            // prepare assembly result
            AssemblyResultSet assemblyResult = new AssemblyResultSet();
            haplotypes.forEach(haplotype -> assemblyResult.add(haplotype));
            Map<String, List<GATKRead>> perSampleReadList = new LinkedHashMap<>();
            List<GATKRead> gtakReads = new LinkedList<>();
            reads.forEach(flowBasedRead -> gtakReads.add(flowBasedRead));
            perSampleReadList.put(sampleNames[0], gtakReads);
            assemblyResult.setPaddedReferenceLoc(loc);
            assemblyResult.setFullReferenceWithPadding(reference.queryAndPrefetch(loc).getBases());

            // computer likelihood
            AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods = likelihoodCalculationEngine.computeReadLikelihoods(
                    assemblyResult, samplesList, perSampleReadList, false);

            // assign
            Map<String, List<GATKRead>> perSampleFilteredReadList = perSampleReadList;
            SimpleInterval genotypingSpan = loc;
            SAMFileHeader readsHeader = readsReader.getHeader();
            Map<Integer,AlleleLikelihoods<GATKRead, Allele>>    genotypeLikelihoods = genotypingEngine.assignGenotypeLikelihoods2(
                    haplotypes,
                    readLikelihoods,
                    perSampleFilteredReadList,
                    assemblyResult.getFullReferenceWithPadding(),
                    assemblyResult.getPaddedReferenceLoc(),
                    genotypingSpan,
                    null,
                    variants,
                    false,
                    hcArgs.maxMnpDistance,
                    readsHeader,
                    false);
            resultWriter.add(loc, genotypeLikelihoods, variants);
        });

        resultWriter.close();
    }
}

