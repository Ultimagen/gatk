package org.ultimagen.featureMapping;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.ultimagen.CommentedTextReader;
import org.ultimagen.FlowTestConstants;

import java.io.File;
import java.io.IOException;

public class FlowFeatureMapperIntegrationTest extends CommandLineProgramTest {

    protected static String    vcTestDir = publicTestDir + FlowTestConstants.FEATURE_MAPPING_DATA_DIR;

    @Test
    public void testBasic() throws IOException {

        final File outputDir = createTempDir("testFlowFeatureMapperTest");
        final File outputFile = new File(outputDir + "/snv_feature_mapper_output.vcf");
        final File expectedFile = new File(vcTestDir + "/snv_feature_mapper_output.vcf");

        final String[] args = new String[] {
                "-R", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "-O", outputFile.getAbsolutePath(),
                "-I", vcTestDir + "/snv_feature_mapper_input.bam",
                "--smith-waterman", "FASTEST_AVAILABLE",
                "--likelihood-calculation-engine", "FlowBased",
                "-mbq", "0",
                "--kmer-size", "10",
                "--copy-attr", "tr",
                "--limit-score", "100",
                "--min-score", "0",
                "--snv-identical-bases", "10",
                "--debug-negatives", "false",
                "--debug-read-name", "150451-BC94-0645901755"
        };

        // run the tool
        runCommandLine(args);  // no assert, just make sure we don't throw

        // make sure we've generated the otuput file
        Assert.assertTrue(outputFile.exists());

        // walk the output and expected files, compare non-comment lines
        try (
                final CommentedTextReader     expectedReader = new CommentedTextReader(expectedFile);
                final CommentedTextReader     outputReader = new CommentedTextReader(outputFile);
        ) {
            String expectedLine;
            String outputLine;
            while ((expectedLine = expectedReader.readLine()) != null) {

                outputLine = outputReader.readLine();
                Assert.assertNotNull(outputLine, "output file contains too few lines");

                Assert.assertEquals(expectedLine, outputLine, "expected and output lines differ");
            }
            outputLine = outputReader.readLine();
            Assert.assertNull(outputLine, "output file contains too many lines");
        }
    }

}
