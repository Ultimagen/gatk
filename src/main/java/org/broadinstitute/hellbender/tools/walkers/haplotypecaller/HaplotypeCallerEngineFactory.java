package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;

public class HaplotypeCallerEngineFactory {

    static public HaplotypeCallerEngine newInstance(final HaplotypeCallerArgumentCollection hcArgs, final AssemblyRegionArgumentCollection assemblyRegionArgs, final boolean createBamOutIndex,
                                                    final boolean createBamOutMD5, final SAMFileHeader readsHeader,
                                                    final ReferenceSequenceFile referenceReader, final VariantAnnotatorEngine annotationEngine) {

        if ( !isRamped(hcArgs) ) {
            return new HaplotypeCallerEngine(hcArgs, assemblyRegionArgs, createBamOutIndex,
                    createBamOutMD5, readsHeader,
                    referenceReader, annotationEngine);
        } else {
            return new ModularHaplotypeCallerEngine(hcArgs, assemblyRegionArgs, createBamOutIndex,
                    createBamOutMD5, readsHeader,
                    referenceReader, annotationEngine);
        }
    }

    private static boolean isRamped(final HaplotypeCallerArgumentCollection hcArgs) {

        return (hcArgs.onRampType != null) || (hcArgs.offRampType != null);
    }
}
