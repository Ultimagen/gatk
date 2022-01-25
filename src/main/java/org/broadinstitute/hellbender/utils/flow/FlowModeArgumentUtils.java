package org.broadinstitute.hellbender.utils.flow;

import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MarkDuplicatesSparkArgumentCollection;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

public class FlowModeArgumentUtils {


    /**
     * the different flow modes, in terms of their parameters and their values
     *
     * NOTE: a parameter value ending with /o is optional - meaning it will not fail the process if it
     * is not existent on the target parameters collection. This allows setting parameters which are
     * specific to only a subset of the tools supporting flow-mode
     */
    public enum FlowModeHC {
        NONE(new String[]{}, null),

        STANDARD(new String[]{
                "mbq", "0",
                "flow-filter-alleles", "true",
                "flow-filter-alleles-sor-threshold", "3",
                "flow-assembly-collapse-hmer-size", "12",
                "flow-matrix-mods", "10,12,11,12",
                "override-fragment-softclip-check", "true",
                "flow-likelihood-parallel-threads", "2",
                "flow-likelihood-optimized-comp", "true",
                "likelihood-calculation-engine", "FlowBased"
        }, null),

        ADVANCED(new String[]{
                "adaptive-pruning/o", "true",
                "pruning-lod-threshold", "3.0",
                "enable-dynamic-read-disqualification-for-genotyping", "true",
                "dynamic-read-disqualification-threshold", "10",
                "apply-frd/o", "true",
                "minimum-mapping-quality", "1",
                "mapping-quality-threshold-for-genotyping/o", "1"
        }, STANDARD);

        final String[] nameValuePairs;
        final FlowModeHC parent;

        FlowModeHC(final String[] nameValuePairs, final FlowModeHC parent) {
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
     * set flow mode defauls for haplotype caller
     */
    public static void setFlowModeHC(CommandLineParser parser, FlowModeHC mode) {
        setArgValues(parser, mode.nameValuePairs);
        if ( mode.parent != null )
            setFlowModeHC(parser, mode.parent);
    }

    /**
     * set flow mode defaults for mark duplicates spark
     */
    public static void setFlowModeMD(CommandLineParser parser) {
        setArgValues(parser, flowModeMD);
    }

    /**
     * common code to set flow mode defauls for anything that uses AssemblyBasedCallerArgumentCollection
     *
     * We made an effort not to override arguments already provided on the commandline
     */
    private static void setArgValues(CommandLineParser parser, String[] nameValuePairs) {

        for ( int i = 0 ; i < nameValuePairs.length ; i += 2 ) {
            if ( !hasBeenSet(parser, cleanParamName(nameValuePairs[i])) ) {
                setValue(parser, nameValuePairs[i], nameValuePairs[i+1]);
            }
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

    private static void setValue(CommandLineParser parser, String alias, String value) {
        if ( parser instanceof CommandLineArgumentParser ) {
            NamedArgumentDefinition namedArg = ((CommandLineArgumentParser)parser).getNamedArgumentDefinitionByAlias(cleanParamName(alias));
            if ( namedArg == null ) {
                if ( isParamOptional(alias) ) {
                    // exit silenently, as it is optional
                    return;
                } else {
                    throw new IllegalArgumentException("alias not found: " + alias);
                }
            }

            PrintStream         ps = new PrintStream(new ByteArrayOutputStream());
            List<String>        values = Arrays.asList(value);
            namedArg.setArgumentValues((CommandLineArgumentParser)parser, ps, values);
        } else {
            throw new IllegalArgumentException("command line parser is not CommandLineArgumentParser");
        }

    }

    private static boolean isParamOptional(final String name) {
        return name.endsWith("/o");
    }

    private static String cleanParamName(final String name) {
        return name.replace("/o", "");
    }

}
