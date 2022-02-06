package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps.AssemblerOffRamp;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps.PostAssemblerOnRamp;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps.PostFilterOnRamp;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps.PreFilterOffRamp;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;

@CommandLineProgramProperties(
        summary = "Call germline SNPs and indels via local re-assembly of haplotypes (ramped version)",
        oneLineSummary = "Call germline SNPs and indels via local re-assembly of haplotypes (ramped version)",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public class RampedHaplotypeCaller extends HaplotypeCaller {

    @ArgumentCollection
    private RampedHaplotypeCallerArgumentCollection rpArgs = new RampedHaplotypeCallerArgumentCollection();

    @Override
    protected HaplotypeCallerEngine buildHaplotypeCallerEngine(final HaplotypeCallerArgumentCollection hcArgs, final AssemblyRegionArgumentCollection assemblyRegionArgs, final boolean createOutputBamIndex, final boolean createOutputBamMD5, final SAMFileHeader headerForReads, final CachingIndexedFastaSequenceFile referenceReader, final VariantAnnotatorEngine variantAnnotatorEngine) {
        return new RampedHaplotypeCallerEngine(hcArgs, assemblyRegionArgs, createOutputBamIndex, createOutputBamMD5, getHeaderForReads(), getReferenceReader(referenceArguments), variantAnnotatorEngine, rpArgs);
    }

}
