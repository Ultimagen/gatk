package org.ultimagen;

import org.testng.Assert;

import java.io.File;
import java.io.IOException;

public class TestFileVerifySame {

    static public void verifySame(final String outputPath, final String expectedPath) throws IOException  {
        verifySame(new File(outputPath), new File(expectedPath));
    }
    static public void verifySame(final File outputPath, final File expectedPath) throws IOException  {

        try (
                final CommentedTextReader expectedReader = new CommentedTextReader(expectedPath);
                final CommentedTextReader outputReader = new CommentedTextReader(outputPath);
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
