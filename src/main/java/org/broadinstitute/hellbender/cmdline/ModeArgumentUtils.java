package org.broadinstitute.hellbender.cmdline;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.*;

public class ModeArgumentUtils {

    protected static final Logger logger = LogManager.getLogger(ModeArgumentUtils.class);

    /**
     * set this mode's defaults - using a set of argValues
     */
    public static void setArgValues(final CommandLineParser parser, final String[] argValues, final String modeName) {
        final Map<String, String>  modifiedArgs = new LinkedHashMap<>();

        for ( int i = 0 ; i < argValues.length ; i += 2 ) {
            if ( !hasBeenSet(parser, argValues[i]) ) {
                if ( setValue(parser, argValues[i], argValues[i+1]) ) {
                    modifiedArgs.put(argValues[i], argValues[i+1]);
                }
            } else {
                logger.info("parameter not set by this mode, as it was already set on the command line: " + argValues[i]);
            }
        }

        logModeNotice(modifiedArgs, modeName);
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
            NamedArgumentDefinition namedArg = ((CommandLineArgumentParser)parser).getNamedArgumentDefinitionByAlias(alias);
            if ( namedArg == null ) {
                throw new IllegalArgumentException("alias not found: " + alias);
            }

            PrintStream         ps = new PrintStream(new ByteArrayOutputStream());
            List<String>        values = Arrays.asList(value);
            namedArg.setArgumentValues((CommandLineArgumentParser)parser, ps, values);
            return true;
        } else {
            throw new IllegalArgumentException("command line parser is not CommandLineArgumentParser");
        }

    }

    private static void logModeNotice(final Map<String, String> modifiedArgs, final String modeName) {
        logger.warn("*************************************************************************");
        logger.warn(String.format("* %-69s *", "--" + modeName + " was enabled"));
        logger.warn("* The following arguments have had their inputs overwritten:            *");
        modifiedArgs.forEach((name, value) -> {
            logger.warn(String.format("* %-69s *", "--" + name + " " + value));
        });
        logger.warn("*                                                                       *");
        logger.warn("* If you would like to run this mode with different inputs for any      *");
        logger.warn("* of the above arguments please manually construct the command or       *");
        logger.warn("* add your specific inputs after the mode argument. This mode           *");
        logger.warn("* will not override inputs explicitly provided.                         *");
        logger.warn("*************************************************************************");
    }
}
