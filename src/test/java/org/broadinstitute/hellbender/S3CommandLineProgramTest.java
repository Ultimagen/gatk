package org.broadinstitute.hellbender;

import org.testng.Assert;
import java.io.IOException;

public class S3CommandLineProgramTest extends CommandLineProgramTest {

    public static final String     S3_BUCKET = "s3://gatk-s3-test";

    @Override
    public String getTestedClassName(){
        return super.getTestedClassName().replaceAll("S3$", "");
    }

    public String getS3TestDir() {
        return S3_BUCKET + "/" + getClass().getSimpleName();
    }

    public String getS3TestFile(String name) {
        return getS3TestDir() + "/" + name;
    }

    public String getGATKS3TestFile(String name) {
        return getS3TestFile(name).replaceAll("^s3://", "s3:///");
    }

    public void execAWS(final String awsCmd)  {
        final String        cmd = "aws " + awsCmd;
        try {
            final Process       process = Runtime.getRuntime().exec(cmd);
            final int           exitCode = process.waitFor();
            Assert.assertTrue(exitCode == 0, "failed to execute: " + cmd);
        } catch (InterruptedException | IOException e ) {
            Assert.fail("interrupted when executing: " + cmd);
        }
    }
}
