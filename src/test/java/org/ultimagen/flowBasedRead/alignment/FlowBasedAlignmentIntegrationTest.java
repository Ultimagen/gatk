package org.ultimagen.flowBasedRead.alignment;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.ultimagen.FlowTestConstants;
import org.ultimagen.TestFileVerifySame;

import java.io.File;
import java.io.IOException;
import java.util.regex.Pattern;

public class FlowBasedAlignmentIntegrationTest extends CommandLineProgramTest {

    private static String    testDir = publicTestDir + "/large";

    @Test
    public void testMatrix() throws IOException {

        final File outputDir = createTempDir("testMatrix");
        final String outputPath = outputDir + "/output.alm";
        final String[] args = new String[] {
                "-R", publicTestDir + "/" + FlowTestConstants.FLOW_BASED_ALIGNMENT_DATA_DIR + "/Homo_sapiens_assembly38.chr9.fasta",
                "-O", outputDir + "/ignored.vcf",
                "-I", testDir + "/input_jukebox_for_test.bam",
                "--smith-waterman", "FASTEST_AVAILABLE",
                "--likelihood-calculation-engine", "FlowBased",
                "-mbq", "0",
                "--kmer-size", "10",
                "--intervals", "chr9:81148694-81177540",
                "--alm-path", outputPath,
                "--alm-interval", "chr9"
        };

        runCommandLine(args);  // no assert, just make sure we don't throw

        // verify that output file has been created
        // walk the output and expected files, compare non-comment lines
        // we are using a specialsed verifier to accomodate for rounding errors
        Assert.assertTrue((new File(outputPath)).exists());
        TestFileVerifySame      tester = new TestFileVerifySame.NearlySameDoubles();
        tester.verifySame(outputPath, publicTestDir + "/" + FlowTestConstants.FLOW_BASED_ALIGNMENT_DATA_DIR + "/input_jukebox_for_test.alm");
    }

    @Override
    public String getTestedToolName() {
        return "HaplotypeCaller";
    }

}
