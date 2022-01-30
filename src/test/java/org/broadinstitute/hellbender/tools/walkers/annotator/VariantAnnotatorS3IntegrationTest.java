package org.broadinstitute.hellbender.tools.walkers.annotator;

import org.broadinstitute.hellbender.S3CommandLineProgramTest;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.File;


public class VariantAnnotatorS3IntegrationTest extends S3CommandLineProgramTest {

    @BeforeTest
    public void setup() throws Exception  {

        // NOTE: this has 'aws' command line dependency!

        // make sure bucket exists
        execAWS("s3 mb " + S3_BUCKET);

        // sync input files
        final String localFolder = getToolTestDataDir();
        execAWS("s3 sync " + localFolder + " " + getS3TestDir());
    }

    @Test
    public void testVCFRead() {

        final String inputVCF = getGATKS3TestFile("input.vcf");
        final File outputVCF = createTempFile("output", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addVCF(inputVCF)
                .addOutput(outputVCF);

        // for now, simply run the tool
        runCommandLine(args.getArgsList());
        Assert.assertTrue(outputVCF.exists(), "failed to create output file: " + outputVCF);
    }

}
