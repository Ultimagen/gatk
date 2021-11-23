package org.broadinstitute.hellbender.exceptions;

public class ToolSuccessfulPrematureExit extends RuntimeException {

    private static final long serialVersionUID = 0L;

    public ToolSuccessfulPrematureExit(final String s) {
        super(s);
    }
}
