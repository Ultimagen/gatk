package org.broadinstitute.hellbender.tools.walkers.groundtruth;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.walkers.variantrecalling.FlowTestConstants;
import org.broadinstitute.hellbender.tools.walkers.variantrecalling.TestFileVerifySame;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class GroundTruthReadsBuilderIntegrationTest extends CommandLineProgramTest {

    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    private static String testDir = publicTestDir + FlowTestConstants.GROUND_TRUTH_DATA_DIR;

    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @Test
    public void testBasic() throws IOException {

        final File outputDir = createTempDir("testGroundTruthTest");
        final File expectedFile = new File(testDir + "/ground_truth_output.csv");
        final File outputFile = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expectedFile : new File(outputDir + "/ground_truth_output.csv");

        final String[] args = new String[] {
                "-R", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "-I", testDir + "/150548_1-UGAv3-2.highconf.q60.chr6_30000000_40000000.cram",
                "--maternal-ref", testDir + "/chr6_HG001_maternal.fa",
                "--paternal-ref", testDir + "/chr6_HG001_paternal.fa",
                "--ancestral-translators-base-path", testDir,
                "--output-csv", outputFile.getAbsolutePath(),
                "--subsampling-ratio", "1.0",
                "--intervals", "chr6:31172223-32980498",
                "--smith-waterman", "FASTEST_AVAILABLE",
                "--likelihood-calculation-engine", "FlowBased",
                "-mbq", "0",
                "--kmer-size", "10",
                "--gt-debug", "false",
                "--output-flow-length", "404",
                "--haplotype-output-padding-size", "0",
                "--fill-trimmed-reads", "false",
                "--fill-softclipped-reads", "false",
                "--max-output-reads", "0",
                "--false-snp-compensation", "true"
        };

        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the otuput file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            (new TestFileVerifySame()).verifySame(outputFile, expectedFile);
        }
    }
}
