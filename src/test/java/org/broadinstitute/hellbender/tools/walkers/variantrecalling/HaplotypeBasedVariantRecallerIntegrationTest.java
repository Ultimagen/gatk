package org.broadinstitute.hellbender.tools.walkers.variantrecalling;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class HaplotypeBasedVariantRecallerIntegrationTest extends CommandLineProgramTest {

    private static String testDir = publicTestDir + FlowTestConstants.VARIANT_CALLING_DATA_DIR;

    @Test
    public void testBasic() throws IOException {

        final File outputDir = createTempDir("testHaplotypeBasedVariantRecallerTest");
        final String outputPath = outputDir + "/output.csv";
        final String[] args = new String[] {
                "--alleles-file-vcf", testDir + "/150292-BC05.vcf.gz",
                "--haplotypes-file-bam", testDir + "/haps_chr5.bam",
                "-I", testDir + "/chr5.bam1.rename.bam",
                "-I", testDir + "/chr5.bam2.rename.bam",
                "--reference", largeFileTestDir + "/Homo_sapiens_assembly38.fasta.gz",
                "--matrix-file-csv", outputPath,
                "--likelihood-calculation-engine", "FlowBased",
                "--phred-scaled-global-read-mismapping-rate", "-1",
                "-L", "chr5:70036625"
        };

        runCommandLine(args);  // no assert, just make sure we don't throw

        // verify that output file has been created
        // walk the output and expected files, compare non-comment lines
        Assert.assertTrue((new File(outputPath)).exists());
        (new TestFileVerifySame.NearlySameDoubles()).verifySame(outputPath, testDir + "/bam1_bam2_output.csv");
    }
}