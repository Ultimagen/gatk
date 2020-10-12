package org.ultimagenomics.variant_calling;

import htsjdk.samtools.*;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.ultimagenomics.flow_based_read.read.FlowBasedRead;
import org.ultimagenomics.haplotype_calling.LHWRefView;

import java.util.*;


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

    private final String    SAMPLE_NAME_DEFAULT = "sm1";

    @ArgumentCollection
    private HaplotypeBasedVariantRecallerArgumentCollection vrArgs = new HaplotypeBasedVariantRecallerArgumentCollection();

    @ArgumentCollection
    private final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();

    public HaplotypeBasedVariantRecaller() {
        super();

        // tool specific default
        hcArgs.likelihoodArgs.phredScaledGlobalReadMismappingRate = -1;
    }

    @Override
    public void traverse() {

        // inits
        final ReadLikelihoodCalculationEngine   likelihoodCalculationEngine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(hcArgs.likelihoodArgs, hcArgs.fbargs);
        final String[]                          sampleNames = {SAMPLE_NAME_DEFAULT};
        final SampleList                        samplesList = new IndexedSampleList(Arrays.asList(sampleNames));
        final HaplotypeCallerGenotypingEngine   genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samplesList, !hcArgs.doNotRunPhysicalPhasing);
        final ReferenceDataSource               reference = ReferenceDataSource.of(referenceArguments.getReferencePath());
        final VariantRecallerResultWriter       resultWriter = new VariantRecallerResultWriter(vrArgs.MATRIX_CSV_FILE);
        final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<VariantContext>(
                vrArgs.ALLELE_VCF_FILE.getAbsolutePath(), null, 0, VariantContext.class);
        final HaplotypeRegionWalker             regionWalker = new HaplotypeRegionWalker(vrArgs, referenceArguments.getReferencePath(), getDefaultCloudPrefetchBufferSize());
        final TrimmedReadsReader                readsReader = new TrimmedReadsReader(vrArgs, referenceArguments.getReferencePath(), getDefaultCloudPrefetchBufferSize());
        final CountingReadFilter                readFilter = makeReadFilter(readsReader.getHeader());
        final SAMSequenceDictionary             samSequenceDictionary = readsReader.getSamSequenceDictionary();
        final List<SimpleInterval>              intervals = hasUserSuppliedIntervals() ? getUserIntervals() : IntervalUtils.getAllIntervalsForReference(samSequenceDictionary);
        readsReader.setReadFilter(readFilter);
        progressMeter.setRecordsBetweenTimeChecks(1);

        // walk regions, as defined by argument
        for ( SimpleInterval region : intervals ) {

            logger.info("region: " + region);
            dataSource.query(region).forEachRemaining(vc -> {

                // walk haplotype (groups) under this variant
                SimpleInterval      vcLoc = new SimpleInterval(vc.getContig(), vc.getStart(), vc.getEnd());
                regionWalker.forBest(vcLoc, bestHaplotypes -> {

                    SimpleInterval  haplotypeSpan = new SimpleInterval(bestHaplotypes.get(0).getGenomeLocation());
                    byte[]          refBases = reference.queryAndPrefetch(haplotypeSpan).getBases();

                    LHWRefView      refView = null;
                    final List<Haplotype> processedHaplotypes = new LinkedList<>();
                    if ( (hcArgs.flowAssemblyCollapseHKerSize > 0)
                                    && LHWRefView.needsCollapsing(refBases, hcArgs.flowAssemblyCollapseHKerSize, logger, false) ) {
                        refView = new LHWRefView(hcArgs.flowAssemblyCollapseHKerSize, refBases, haplotypeSpan, logger, false);
                        processedHaplotypes.addAll(refView.uncollapseHaplotypesByRef(bestHaplotypes, false, true, refBases));
                    }
                    else
                        processedHaplotypes.addAll(bestHaplotypes);

                    // get reads overlapping haplotypes
                    Collection<FlowBasedRead> reads = readsReader.getReads(haplotypeSpan, vcLoc);
                    List<VariantContext>      variants = new LinkedList<>(Arrays.asList(vc));
                    if ( logger.isDebugEnabled() )
                        logger.debug(String.format("vcLoc %s, haplotypeSpan: %s, %d haplotypes, %d reads",
                                vcLoc.toString(), haplotypeSpan.toString(), processedHaplotypes.size(), reads.size(), variants.size()));
                    progressMeter.update(vcLoc);

                    // prepare assembly result
                    AssemblyResultSet assemblyResult = new AssemblyResultSet();
                    processedHaplotypes.forEach(haplotype -> assemblyResult.add(haplotype));
                    Map<String, List<GATKRead>> perSampleReadList = new LinkedHashMap<>();
                    List<GATKRead> gtakReads = new LinkedList<>();
                    reads.forEach(flowBasedRead -> gtakReads.add(flowBasedRead));
                    perSampleReadList.put(sampleNames[0], gtakReads);
                    AssemblyRegion regionForGenotyping = new AssemblyRegion(haplotypeSpan, 0, readsReader.getHeader());
                    assemblyResult.setPaddedReferenceLoc(haplotypeSpan);
                    assemblyResult.setFullReferenceWithPadding(refBases);
                    assemblyResult.setRegionForGenotyping(regionForGenotyping);


                    // computer likelihood
                    AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods = likelihoodCalculationEngine.computeReadLikelihoods(
                            assemblyResult, samplesList, perSampleReadList, false);

                    // assign
                    Map<String, List<GATKRead>> perSampleFilteredReadList = perSampleReadList;
                    SAMFileHeader readsHeader = readsReader.getHeader();
                    Map<Integer, AlleleLikelihoods<GATKRead, Allele>> genotypeLikelihoods = genotypingEngine.assignGenotypeLikelihoods2(
                            processedHaplotypes,
                            readLikelihoods,
                            perSampleFilteredReadList,
                            assemblyResult.getFullReferenceWithPadding(),
                            assemblyResult.getPaddedReferenceLoc(),
                            regionForGenotyping.getSpan(),
                            null,
                            variants,
                            false,
                            hcArgs.maxMnpDistance,
                            readsHeader,
                            false);
                    resultWriter.add(haplotypeSpan, genotypeLikelihoods, variants, assemblyResult);
                });
            });
        }

        resultWriter.close();
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return HaplotypeCallerEngine.makeStandardHCReadFilters();
    }

    public CountingReadFilter makeReadFilter(SAMFileHeader samFileHeader){
        final GATKReadFilterPluginDescriptor readFilterPlugin =
                getCommandLineParser().getPluginDescriptor(GATKReadFilterPluginDescriptor.class);
        return  readFilterPlugin.getMergedCountingReadFilter(samFileHeader);
    }

}

