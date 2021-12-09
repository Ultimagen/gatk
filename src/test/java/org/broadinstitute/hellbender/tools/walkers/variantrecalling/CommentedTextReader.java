package org.broadinstitute.hellbender.tools.walkers.variantrecalling;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class CommentedTextReader implements AutoCloseable {
    final boolean         ignoreComments;
    final String          commentPrefix;
    final BufferedReader reader;

    public CommentedTextReader(final File file, final boolean ignoreComments, final String commentPrefix) throws IOException {
        this.ignoreComments = ignoreComments;
        this.commentPrefix = commentPrefix;
        reader = new BufferedReader(new FileReader(file));
    }
    public CommentedTextReader(final File file) throws IOException {
        this(file, true, "#");
    }

    public String readLine() throws IOException {
        String      line;
        while ( (line = reader.readLine()) != null ) {
            if ( !ignoreComments || !line.startsWith(commentPrefix) ) {
                return line;
            }
        }
        return null;
    }

    @Override
    public void close() throws IOException {
        reader.close();
    }
}
