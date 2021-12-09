package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps.*;

import java.io.IOException;
import java.util.*;

/**
 * The core engine for the HaplotypeCaller that does all of the actual work of the tool.
 *
 * Usage:
 * -Pass the HaplotypeCaller args into the constructor, which will initialize the HC engine completely.
 * -Get the appropriate VCF or GVCF writer (depending on our arguments) from {@link #makeVCFWriter}
 * -Write the appropriate VCF header via {@link #writeHeader}
 * -Repeatedly call {@link #isActive} to identify active vs. inactive regions
 * -Repeatedly call {@link #callRegion} to call variants in each region, and add them to your writer
 * -When done, call {@link #shutdown}. Close the writer you got from {@link #makeVCFWriter} yourself.
 */
public abstract class RampedHaplotypeCallerEngine extends HaplotypeCallerEngine {

    private static final Logger logger = LogManager.getLogger(RampedHaplotypeCallerEngine.class);

    // on-off ramp, if present
    protected PreFilterOffRamp preFilterOffRamp = null;
    protected PostFilterOnRamp postFilterOnRamp = null;
    protected AssemblerOffRamp preAssemblerOffRamp = null;
    protected AssemblerOffRamp postAssemblerOffRamp = null;
    protected PostAssemblerOnRamp postAssemblerOnRamp = null;

    public RampedHaplotypeCallerEngine(final HaplotypeCallerArgumentCollection hcArgs, AssemblyRegionArgumentCollection assemblyRegionArgs, boolean createBamOutIndex,
                                       boolean createBamOutMD5, final SAMFileHeader readsHeader,
                                       ReferenceSequenceFile referenceReader, VariantAnnotatorEngine annotationEngine) {

        super(hcArgs, assemblyRegionArgs, createBamOutIndex,
                createBamOutMD5, readsHeader,
                referenceReader, annotationEngine);

        buildRamps();
    }

    /**
     * Shutdown this HC engine, closing resources as appropriate
     */
    public void shutdown() {
        super.shutdown();
        tearRamps();
    }

    public void buildRamps() {
        try {

            if ( hcArgs.offRampType != null && hcArgs.offRampType != AssemblyBasedCallerArgumentCollection.OffRampTypeEnum.NONE) {
                if ( hcArgs.offRampFile == null )
                    throw new RuntimeException("rampFile must be specified");

                // create ramp
                switch ( hcArgs.offRampType ) {
                    case NONE:
                        break;
                    case PRE_FILTER_OFF:
                        preFilterOffRamp = new PreFilterOffRamp(hcArgs.offRampFile);
                        break;
                    case POST_ASSEMBLER_OFF:
                        postAssemblerOffRamp = new AssemblerOffRamp(hcArgs.offRampFile);
                        break;
                    case PRE_ASSEMBLER_OFF:
                        preAssemblerOffRamp = new AssemblerOffRamp(hcArgs.offRampFile);
                        break;
                }
            }

            if ( hcArgs.onRampType != null && hcArgs.onRampType != AssemblyBasedCallerArgumentCollection.OnRampTypeEnum.NONE ) {
                if ( hcArgs.onRampFile == null )
                    throw new RuntimeException("rampFile must be specified");

                // create ramp
                switch ( hcArgs.onRampType ) {
                    case NONE:
                        break;
                    case POST_FILTER_ON:
                        postFilterOnRamp = new PostFilterOnRamp(hcArgs.onRampFile);
                        break;
                    case POST_ASSEMBLER_ON:
                        postAssemblerOnRamp = new PostAssemblerOnRamp(hcArgs.onRampFile);
                        break;
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private void tearRamps() {
        try {
            for (RampBase ramp : Arrays.asList(preFilterOffRamp, postFilterOnRamp,
                    preAssemblerOffRamp, postAssemblerOffRamp, postAssemblerOnRamp) ) {
                if ( ramp != null ) {
                    ramp.close();
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
