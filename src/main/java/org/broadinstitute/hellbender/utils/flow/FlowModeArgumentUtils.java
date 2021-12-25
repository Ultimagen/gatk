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

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

public class FlowModeArgumentUtils {

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
                "adaptive-pruning", "true",
                "pruning-lod-threshold", "3.0",
                "enable-dynamic-read-disqualification-for-genotyping", "true",
                "dynamic-read-disqualification-threshold", "10",
                "apply-frd", "true",
                "minimum-mapping-quality", "1",
                "mapping-quality-threshold-for-genotyping", "1"
        }, STANDARD);

        final String[] nameValuePairs;
        final FlowModeHC parent;

        FlowModeHC(final String[] nameValuePairs, final FlowModeHC parent) {
            this.nameValuePairs = nameValuePairs;
            this.parent = parent;
        }
    } ;

    private static final String[] flowModeMD = {
            MarkDuplicatesSparkArgumentCollection.FLOW_END_LOCATION_SIGNIFICANT_LONG_NAME, "true",
            MarkDuplicatesSparkArgumentCollection.FLOW_USE_CLIPPED_LOCATIONS_LONG_NAME, "true",
            MarkDuplicatesSparkArgumentCollection.ENDS_READ_UNCERTAINTY_LONG_NAME, "true",
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
            if ( !hasBeenSet(parser, nameValuePairs[i]) ) {
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
            NamedArgumentDefinition namedArg = ((CommandLineArgumentParser)parser).getNamedArgumentDefinitionByAlias(alias);
            if ( namedArg == null ) {
                throw new IllegalArgumentException("alias not found: " + alias);
            }

            PrintStream         ps = new PrintStream(new ByteArrayOutputStream());
            List<String>        values = Arrays.asList(value);
            namedArg.setArgumentValues((CommandLineArgumentParser)parser, ps, values);
        } else {
            throw new IllegalArgumentException("command line parser is not CommandLineArgumentParser");
        }

    }

}
