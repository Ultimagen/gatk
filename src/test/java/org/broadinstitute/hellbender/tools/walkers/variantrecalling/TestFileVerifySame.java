package org.broadinstitute.hellbender.tools.walkers.variantrecalling;

import org.testng.Assert;

import java.io.File;
import java.io.IOException;
import java.util.regex.Pattern;

public class TestFileVerifySame {

    public static class NearlySameDoubles extends TestFileVerifySame {

        private static Pattern   isDouble = Pattern.compile("[-+]?[0-9]*\\.?[0-9]+");
        private double           maxDoubleDeltaRatioAllowed = 0.0001;

        @Override
        public void assertEqual(String expectedLine, String outputLine) {

            // if same, escape quitely
            if ( expectedLine.equals(outputLine) )
                return;

            // break lines into tokens that must be the same or close enough (as doubles)
            String[]       expectedToks = expectedLine.replace('\t', ' ').split(" ");
            String[]       outputToks = outputLine.replace('\t', ' ').split(" ");
            if ( expectedToks.length != outputToks.length ) {
                super.assertEqual(expectedLine, outputLine);
            }
            for ( int i = 0 ; i < expectedToks.length ; i++ ) {

                // same?
                if ( expectedToks[i].equals(outputToks[i]) )
                    continue;

                // must be numbers
                if ( !isDouble.matcher(expectedToks[i]).matches()
                        || !isDouble.matcher(outputToks[i]).matches() ) {
                    // make fail
                    super.assertEqual(expectedLine, outputLine);
                }

                double      expected = Double.parseDouble(expectedToks[i]);
                double      output = Double.parseDouble(outputToks[i]);

                // if exactly the same, its ok of course
                if ( expected == output ) {
                    continue;
                }

                // we don't want zeros where they are not expected
                if ( (expected == 0.0) ^ (output == 0.0) ) {
                    // fail - one is zero while the other not
                    Assert.assertEquals(expected, output);
                }

                // we don't want them to be on another side of the zero as well
                if ( (expected > 0.0) ^ (output > 0.0) ) {
                    // fail - one is zero while the other not
                    Assert.assertEquals(expected, output);
                }

                // compute delta ratio, that that below threshold
                double      deltaRatio = Math.abs((expected - output) / ((expected + output) / 2));
                Assert.assertTrue(deltaRatio < maxDoubleDeltaRatioAllowed,
                        String.format("deltaRatio exceeded between %s and %s", expectedToks[i], outputToks[i]));
            }
        }
    }

    public void verifySame(final String outputPath, final String expectedPath) throws IOException  {
        verifySame(new File(outputPath), new File(expectedPath));
    }
    public void verifySame(final File outputPath, final File expectedPath) throws IOException  {

        try (
                final CommentedTextReader expectedReader = new CommentedTextReader(expectedPath);
                final CommentedTextReader outputReader = new CommentedTextReader(outputPath);
        ) {
            String expectedLine;
            String outputLine;
            while ((expectedLine = expectedReader.readLine()) != null) {

                outputLine = outputReader.readLine();
                Assert.assertNotNull(outputLine, "output file contains too few lines");

                assertEqual(expectedLine, outputLine);
            }
            outputLine = outputReader.readLine();
            Assert.assertNull(outputLine, "output file contains too many lines");
        }
    }

    public void assertEqual(String expectedLine, String outputLine) {
        Assert.assertEquals(expectedLine, outputLine, "expected and output lines differ");
    }
};

