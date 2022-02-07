package org.broadinstitute.hellbender;

import com.amazonaws.regions.Regions;
import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.testng.Assert;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

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
            logger.info("using bucket (env): " + bucketName);
        } else {
            bucketName = S3_BUCKET_PREFIX + s3.getS3AccountOwner().getDisplayName().replaceAll("\\W", "");
            logger.info("using bucket (auto): " + bucketName);

            // make sure bucket exists
            execAws(Arrays.asList("aws", "s3", "mb", getBucketName()));
        }

    }

    public String getS3ToolTestDataDir(){
        return "src/test/resources/large/s3/" + getTestedClassName() + "/";
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
        execAws(Arrays.asList("aws", "s3", "sync", getToolTestDataDir(), getS3TestDir()));
    }

    public void execAws(final List<String> awsCmd)  {
        logger.info("execAws: " + awsCmd);
        try {
            final Process       process = Runtime.getRuntime().exec(awsCmd.toArray(new String[0]));
            final int           exitCode = process.waitFor();
            Assert.assertTrue(exitCode == 0, "failed to execute: " + awsCmd);
        } catch (InterruptedException | IOException e ) {
            Assert.fail("interrupted when executing: " + awsCmd);
        }
    }

    public String getBucketName() {
        return bucketName;
    }
}
