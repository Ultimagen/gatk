package org.ultimagen.variantRecalling;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.ultimagen.CommentedTextReader;
import org.ultimagen.FlowTestConstants;

import java.io.File;
import java.io.IOException;

public class HaplotypeBasedVariantRecallerIntegrationTest extends CommandLineProgramTest {

    protected static String    vcTestDir = publicTestDir + FlowTestConstants.VARIANT_CALLING_DATA_DIR;

    @Test
    public void testBasic() throws IOException {

        final File outputDir = createTempDir("testHaplotypeBasedVariantRecallerTest");
        final String outputPath = outputDir + "/output.csv";
        final String[] args = new String[] {
                "--alleles-file-vcf", vcTestDir + "/150292-BC05.vcf.gz",
                "--haplotypes-file-bam", vcTestDir + "/haps_chr5.bam",
                "-I", vcTestDir + "/chr5.bam1.rename.bam",
                "-I", vcTestDir + "/chr5.bam2.rename.bam",
                "--reference", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "--matrix-file-csv", outputPath,
                "--likelihood-calculation-engine", "FlowBased",
                "--phred-scaled-global-read-mismapping-rate", "-1",
                "-L", "chr5:70036625"
        };

        runCommandLine(args);  // no assert, just make sure we don't throw

        // verify that output file has been created
        Assert.assertTrue((new File(outputPath)).exists());

        // walk the output and expected files, compare non-comment lines
        final String expectedPath = outputDir + "/bam1_bam2_output.csv";
        try (
                final CommentedTextReader expectedReader = new CommentedTextReader(new File(outputPath));
                final CommentedTextReader outputReader = new CommentedTextReader(new File(outputPath));
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
