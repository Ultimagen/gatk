package org.broadinstitute.hellbender.utils.flow;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;
import org.broadinstitute.hellbender.tools.FlowBasedAlignmentArgumentCollection;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.*;

public class FlowModeArgumentUtils {

    protected static final Logger logger = LogManager.getLogger(FlowModeArgumentUtils.class);

    private static final String OPTIONAL_SUFFIX = "/o";

    /**
     * the different flow modes, in terms of their parameters and their values
     *
     * NOTE: a parameter value ending with /o is optional - meaning it will not fail the process if it
     * is not existent on the target parameters collection. This allows setting parameters which are
     * specific to only a subset of the tools supporting flow-mode
     */
    public enum FlowMode {
        NONE(new String[]{}, null),

        STANDARD(new String[]{
                AssemblyBasedCallerArgumentCollection.MIN_BASE_QUALITY_SCORE_SHORT_NAME, "0",
                AssemblyBasedCallerArgumentCollection.FILTER_ALLELES, "true",
                AssemblyBasedCallerArgumentCollection.FILTER_ALLELES_SOR_THRESHOLD, "3",
                AssemblyBasedCallerArgumentCollection.FLOW_ASSEMBLY_COLLAPSE_HMER_SIZE_LONG_NAME, "12",
                FlowBasedArgumentCollection.FLOW_MATRIX_MODS_LONG_NAME, "10,12,11,12",
                AssemblyBasedCallerArgumentCollection.OVERRIDE_FRAGMENT_SOFTCLIP_CHECK_LONG_NAME, "true",
                FlowBasedAlignmentArgumentCollection.FLOW_LIKELIHOOD_PARALLEL_THREADS_LONG_NAME, "2",
                FlowBasedAlignmentArgumentCollection.FLOW_LIKELIHOOD_OPTIMIZED_COMP_LONG_NAME, "true",
                LikelihoodEngineArgumentCollection.LIKELIHOOD_CALCULATION_ENGINE_FULL_NAME, "FlowBased"
        }, null),

        ADVANCED(new String[]{
                HaplotypeCallerReadThreadingAssemblerArgumentCollection.ADAPTIVE_PRUNING_LONG_NAME + OPTIONAL_SUFFIX, "true",
                ReadThreadingAssemblerArgumentCollection.PRUNING_LOD_THRESHOLD_LONG_NAME, "3.0",
                LikelihoodEngineArgumentCollection.ENABLE_DYNAMIC_READ_DISQUALIFICATION_FOR_GENOTYPING_LONG_NAME, "true",
                LikelihoodEngineArgumentCollection.DYNAMIC_READ_DISQUALIFICATION_THRESHOLD_LONG_NAME, "10",
                HaplotypeCallerArgumentCollection.APPLY_FRD_LONG_NAME + OPTIONAL_SUFFIX, "true",
                ReadFilterArgumentDefinitions.MINIMUM_MAPPING_QUALITY_NAME, "1",
                HaplotypeCallerArgumentCollection.MAPPING_QUALITY_THRESHOLD_FOR_GENOTYPING_LONG_NAME + OPTIONAL_SUFFIX, "1"
        }, STANDARD);

        final String[] nameValuePairs;
        final FlowMode parent;

        FlowMode(final String[] nameValuePairs, final FlowMode parent) {
            this.nameValuePairs = nameValuePairs;
            this.parent = parent;
        }
    } ;

    private static final String[] flowModeMD = {
            MarkDuplicatesSparkArgumentCollection.SINGLE_END_READS_END_POSITION_SIGNIFICANT, "true",
            MarkDuplicatesSparkArgumentCollection.SINGLE_END_READS_CLIPPING_IS_END_LONG_NAME, "true",
            MarkDuplicatesSparkArgumentCollection.FLOW_END_POS_UNCERTAINTY_LONG_NAME, "1",
            MarkDuplicatesSparkArgumentCollection.FLOW_SKIP_START_HOMOPOLYMERS_LONG_NAME, "0"
    };


    /**
     * set flow mode defauls for haplotype caller. returns map of args actually modified
     */
    public static Map<String, String> setFlowMode(final CommandLineParser parser, final FlowMode mode) {
        final Map<String, String>  modifiedArgs = setArgValues(parser, mode);
        logFlowModeNotice(modifiedArgs);
        return modifiedArgs;
    }

    /**
     * set flow mode defaults for mark duplicates spark
     */
    public static Map<String, String> setFlowModeMD(final CommandLineParser parser) {
        final Map<String, String>  modifiedArgs = setArgValues(parser, flowModeMD);
        logFlowModeNotice(modifiedArgs);
        return modifiedArgs;
    }

    /**
     * common code to set flow mode defauls for anything that uses AssemblyBasedCallerArgumentCollection
     *
     * We made an effort not to override arguments already provided on the commandline
     */
    private static Map<String, String> setArgValues(final CommandLineParser parser, final String[] nameValuePairs) {

        final Map<String, String>  modifiedArgs = new LinkedHashMap<>();

        for ( int i = 0 ; i < nameValuePairs.length ; i += 2 ) {
            if ( !hasBeenSet(parser, cleanParamName(nameValuePairs[i])) ) {
                if ( setValue(parser, nameValuePairs[i], nameValuePairs[i+1]) ) {
                    modifiedArgs.put(cleanParamName(nameValuePairs[i]), nameValuePairs[i+1]);
                }
            }
        }

        return modifiedArgs;
    }

    private static Map<String, String> setArgValues(final CommandLineParser parser, final FlowMode mode) {

        final Map<String, String>  modifiedArgs =setArgValues(parser, mode.nameValuePairs);
        if ( mode.parent != null ) {
            modifiedArgs.putAll(setArgValues(parser, mode.parent));
        }

        return modifiedArgs;
    }


    private static boolean hasBeenSet(final CommandLineParser parser, final String alias) {

        if ( parser instanceof CommandLineArgumentParser ) {
            NamedArgumentDefinition namedArg = ((CommandLineArgumentParser)parser).getNamedArgumentDefinitionByAlias(alias);

            return (namedArg != null) ? namedArg.getHasBeenSet() : false;
        } else {
            return false;
        }
    }

    private static boolean setValue(final CommandLineParser parser, final String alias, final String value) {
        if ( parser instanceof CommandLineArgumentParser ) {
            NamedArgumentDefinition namedArg = ((CommandLineArgumentParser)parser).getNamedArgumentDefinitionByAlias(cleanParamName(alias));
            if ( namedArg == null ) {
                if ( isParamOptional(alias) ) {
                    // exit silenently, as it is optional
                    return false;
                } else {
                    throw new IllegalArgumentException("alias not found: " + alias);
                }
            }

            PrintStream         ps = new PrintStream(new ByteArrayOutputStream());
            List<String>        values = Arrays.asList(value);
            namedArg.setArgumentValues((CommandLineArgumentParser)parser, ps, values);
            return true;
        } else {
            throw new IllegalArgumentException("command line parser is not CommandLineArgumentParser");
        }

    }

    private static boolean isParamOptional(final String name) {
        return name.endsWith(OPTIONAL_SUFFIX);
    }

    private static String cleanParamName(final String name) {
        return name.replace(OPTIONAL_SUFFIX, "");
    }

    private static void logFlowModeNotice(final Map<String, String> modifiedArgs) {
        logger.warn("*************************************************************************");
        logger.warn("* --flow-mode was enabled                                               *");
        logger.warn("* The following arguments have had their inputs overwritten:            *");
        modifiedArgs.forEach((name, value) -> {
            logger.warn(String.format("* %-69s *", "--" + name + " " + value));
        });
        logger.warn("*                                                                       *");
        logger.warn("* If you would like to run flow mode with different inputs for any      *");
        logger.warn("* of the above arguments please manually construct the command or       *");
        logger.warn("* add your specific inputs after the flow-mode argument. Flow mode      *");
        logger.warn("* will not override inputs explicitly provided.                         *");
        logger.warn("*************************************************************************");
    }
}
