package org.ultimagenomics.variant_calling;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class HaplotypeBasedVariantRecallerTest extends CommandLineProgramTest {

    protected static String    vcTestDir = publicTestDir + Const.DATA_DIR;
    protected static String    refDir = gatkDirectory + "testdata/ref";

    @Test
    public void testBasic() throws IOException {

        final File outputDir = createTempDir("testHaplotypeBasedVariantRecallerTest");
        final String[] args = new String[] {
                "--alleles-file-vcf", vcTestDir + "/150292-BC05.vcf.gz",
                "--haplotypes-file-bam", vcTestDir + "/haps_chr5.bam",
                "-I", vcTestDir + "/chr5.bam1.rename.bam",
                "-I", vcTestDir + "/chr5.bam2.rename.bam",
                "--reference", refDir + "/Homo_sapiens_assembly38.fasta",
                "--matrix-file-csv", outputDir + "/output.csv",
                "--likelihood-calculation-engine", "FlowBased",
                "--phred-scaled-global-read-mismapping-rate", "-1",
                "-L", "chr5:70036625"
        };

        runCommandLine(args);  // no assert, just make sure we don't throw
    }
}
