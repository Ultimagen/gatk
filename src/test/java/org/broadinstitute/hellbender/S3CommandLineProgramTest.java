package org.broadinstitute.hellbender;

import com.amazonaws.regions.Regions;
import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.testng.Assert;
import java.io.IOException;

public class S3CommandLineProgramTest extends CommandLineProgramTest {

    private static final Logger logger = LogManager.getLogger(S3CommandLineProgramTest.class);

    private static final String     S3_BUCKET_PREFIX = "s3://gatk-s3-test-";
    private static final String     S3_BUCKET_ENV_NAME = "GATK_S3_TEST_BUCKET";

    private final String            bucketName;

    public S3CommandLineProgramTest() {
        super();

        // find out current aws user id
        AmazonS3 s3 = AmazonS3ClientBuilder.standard().withRegion(Regions.DEFAULT_REGION).build();

        // estabish bucket name
        if ( System.getenv(S3_BUCKET_ENV_NAME) != null ) {
            bucketName = System.getenv(S3_BUCKET_ENV_NAME);
        } else {
            bucketName = S3_BUCKET_PREFIX + s3.getS3AccountOwner().getDisplayName().replaceAll("\\W", "");
        }
        logger.info("using bucket: " + bucketName);

        // make sure bucket exists
        execAws("s3 mb " + getBucketName());
    }

    @Override
    public String getTestedClassName(){
        return super.getTestedClassName().replaceAll("S3$", "");
    }

    public String getS3TestDir() {
        return getBucketName() + "/" + getClass().getSimpleName();
    }

    public String getS3TestFile(String name) {
        return getS3TestDir() + "/" + name;
    }

    public String getGATKS3TestFile(String name) {
        return getS3TestFile(name).replaceAll("^s3://", "s3:///");
    }

    public void syncToolTestData() {
        execAws("s3 sync " + getToolTestDataDir() + " " + getS3TestDir());
    }

    public void execAws(final String awsCmd)  {
        final String        cmd = "aws " + awsCmd;
        try {
            final Process       process = Runtime.getRuntime().exec(cmd);
            final int           exitCode = process.waitFor();
            Assert.assertTrue(exitCode == 0, "failed to execute: " + cmd);
        } catch (InterruptedException | IOException e ) {
            Assert.fail("interrupted when executing: " + cmd);
        }
    }

    public String getBucketName() {
        return bucketName;
    }
}
