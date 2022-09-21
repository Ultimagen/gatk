package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.tools.walkers.variantrecalling.FlowTestConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class SelectInformativeReadsIntegrationTest extends CommandLineProgramTest {

    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    private static String testDir = publicTestDir + FlowTestConstants.FEATURE_MAPPING_DATA_DIR;

    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @Test
    public void testBasic() throws IOException {

        final File outputDir = createTempDir("test" + getClass().getSimpleName());
        final String outputName = "/" + getClass().getSimpleName() + "_output.sam";
        final File expectedFile = new File(testDir + outputName);
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + outputName);
        final File input = new File(largeFileTestDir, "input_jukebox_for_test.bam");
        final File alleleFile = new File(largeFileTestDir, "input_jukebox_for_test_dbSNP.vcf");

        final String[] args = new String[] {
                "-I", input.getAbsolutePath(),
                "-O", outputFile.getAbsolutePath(),
                "--intervals", "chr9:81149486-81177047",
                "--allele-file", alleleFile.getAbsolutePath(),
                "--min-ref-allele-distance", "0.2",
                "--max-abs-ref-score", "0.5"
        };

        // run the tool
        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the otuput file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(outputFile, expectedFile, "@");
        }
    }

}
