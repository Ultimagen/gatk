package org.broadinstitute.hellbender.utils.flow;

import com.google.common.collect.Lists;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

public class FlowModeArgumentUtils {

    /**
     * set flow mode defauls for haplotype caller
     */
    public static void setModeDefaults(HaplotypeCallerArgumentCollection args) {

        set(args);
        args.fbargs.flowLikelihoodOptimizedComp = true;
    }

    /**
     * set flow mode defaults for mutect2
     */
    public static void setModeDefaults(M2ArgumentCollection args) {
        set(args);
        args.fbargs.flowLikelihoodOptimizedComp = true;
    }

    /**
     * common code to set flow mode defauls for anything that uses AssemblyBasedCallerArgumentCollection
     */
    private static void set(AssemblyBasedCallerArgumentCollection args) {

        args.smithWatermanImplementation = SmithWatermanAligner.Implementation.FASTEST_AVAILABLE;
        args.likelihoodArgs.likelihoodEngineImplementation = ReadLikelihoodCalculationEngine.Implementation.FlowBased;
        args.minBaseQualityScore = 0;
        args.assemblerArgs.kmerSizes = Lists.newArrayList(10);
        args.flowMatrixMods = "10,12,11,12";
        args.flowAssemblyCollapseHKerSize = 12;
        args.prefilterSorThreshold = 40;
        args.filterAlleles = true;
        args.filterLoneAlleles = true;
    }

    /**
     * set flow mode defaults for mark duplicates spark
     */
    public static void setModeDefaults(MarkDuplicatesSparkArgumentCollection args) {

        /**
         * TODO: fill these when values become known (or delete if no default needs to be set)
         *
         * NOTE: current values are the instance defaults
         */
        args.FLOW_QUALITY_SUM_STRATEGY = false;
        args.FLOW_END_LOCATION_SIGNIFICANT = false;
        args.ENDS_READ_UNCERTAINTY = 0;
        args.FLOW_USE_CLIPPED_LOCATIONS = false;
        args.FLOW_SKIP_START_HOMOPOLYMERS = 0;
        args.FLOW_Q_IS_KNOWN_END = false;
    }
}
