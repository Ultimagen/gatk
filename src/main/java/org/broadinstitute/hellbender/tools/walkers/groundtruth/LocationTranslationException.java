package org.broadinstitute.hellbender.tools.walkers.groundtruth;

public class LocationTranslationException extends RuntimeException {

    static final private long        serialVersionUID = 0;

    LocationTranslationException(String msg) {

        super(msg);
    }
}
