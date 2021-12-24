package org.broadinstitute.hellbender.utils.flow;

import com.google.common.collect.Lists;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.FlowBasedAlignmentArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.LikelihoodEngineArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

public class FlowModeArgumentUtils {

    /**
     * set flow mode defauls for haplotype caller
     */
    public static void setModeDefaults(CommandLineParser parser, HaplotypeCallerArgumentCollection args) {

        set(parser, args);
        if ( !hasBeenSet(parser, FlowBasedAlignmentArgumentCollection.FLOW_LIKELIHOOD_OPTIMIZED_COMP) ) {
            args.fbargs.flowLikelihoodOptimizedComp = true;
        }
    }

    /**
     * set flow mode defaults for mutect2
     */
    public static void setModeDefaults(CommandLineParser parser, M2ArgumentCollection args) {
        set(parser, args);
        if ( !hasBeenSet(parser, FlowBasedAlignmentArgumentCollection.FLOW_LIKELIHOOD_OPTIMIZED_COMP) ) {
            args.fbargs.flowLikelihoodOptimizedComp = true;
        }
    }

    /**
     * common code to set flow mode defauls for anything that uses AssemblyBasedCallerArgumentCollection
     *
     * We made an effort not to override arguments already provided on the commandline
     */
    private static void set(CommandLineParser parser, AssemblyBasedCallerArgumentCollection args) {

        if ( !hasBeenSet(parser, AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_LONG_NAME) ) {
            args.smithWatermanImplementation = SmithWatermanAligner.Implementation.JAVA;
        }

        if ( !hasBeenSet(parser, LikelihoodEngineArgumentCollection.LIKELIHOOD_CALCULATION_ENGINE_FULL_NAME) ) {
            args.likelihoodArgs.likelihoodEngineImplementation = ReadLikelihoodCalculationEngine.Implementation.FlowBased;
        }

        if ( !hasBeenSet(parser, AssemblyBasedCallerArgumentCollection.MIN_BASE_QUALITY_SCORE_LONG_NAME) ) {
            args.minBaseQualityScore = 0;
        }

        if ( !hasBeenSet(parser, AssemblyBasedCallerArgumentCollection.FLOW_MATRIX_MODS_LONG_NAME) ) {
            args.flowMatrixMods = "10,12,11,12";
        }

        if ( !hasBeenSet(parser, AssemblyBasedCallerArgumentCollection.FLOW_ASSEMBLY_COLLAPSE_HMER_SIZE_LONG_NAME) ) {
            args.flowAssemblyCollapseHKerSize = 12;
        }

        if ( !hasBeenSet(parser, AssemblyBasedCallerArgumentCollection.FILTER_ALLELES) ) {
            args.filterAlleles = true;
        }

        if ( !hasBeenSet(parser, AssemblyBasedCallerArgumentCollection.FILTER_ALLELES_FILTER_LONE_ALLELES) ) {
            args.filterLoneAlleles = false;
        }

        if ( !hasBeenSet(parser, AssemblyBasedCallerArgumentCollection.FILTER_ALLELES_SOR_THRESHOLD) ) {
            args.prefilterSorThreshold = 3;
        }
    }

    /**
     * set flow mode defaults for mark duplicates spark
     */
    public static void setModeDefaults(CommandLineParser parser, MarkDuplicatesSparkArgumentCollection args) {

        /**
         * TODO: fill these when values become known (or delete if no default needs to be set)
         *
         * NOTE: current values are the instance defaults
         */
        if ( !hasBeenSet(parser, MarkDuplicatesSparkArgumentCollection.FLOW_END_LOCATION_SIGNIFICANT_LONG_NAME) ) {
            args.FLOW_END_LOCATION_SIGNIFICANT = true;
        }

        if ( !hasBeenSet(parser, MarkDuplicatesSparkArgumentCollection.FLOW_USE_CLIPPED_LOCATIONS_LONG_NAME) ) {
            args.FLOW_USE_CLIPPED_LOCATIONS = true;
        }

        if ( !hasBeenSet(parser, MarkDuplicatesSparkArgumentCollection.ENDS_READ_UNCERTAINTY_LONG_NAME) ) {
            args.ENDS_READ_UNCERTAINTY = 1;
        }

        if ( !hasBeenSet(parser, MarkDuplicatesSparkArgumentCollection.FLOW_SKIP_START_HOMOPOLYMERS_LONG_NAME) ) {
            args.FLOW_SKIP_START_HOMOPOLYMERS = 0;
        }
    }

    private static boolean hasBeenSet(CommandLineParser parser, String alias) {

        if ( parser instanceof CommandLineArgumentParser ) {
            NamedArgumentDefinition namedArg = ((CommandLineArgumentParser)parser).getNamedArgumentDefinitionByAlias(alias);

            return (namedArg != null) ? namedArg.getHasBeenSet() : false;
        } else {
            return false;
        }

    }


}
