package org.ultimagenomics.ramps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import htsjdk.samtools.util.Locatable;
import org.json.JSONObject;

import java.io.File;
import java.io.IOException;

public abstract class RampBase {

    protected static final Logger logger = LogManager.getLogger(RampBase.class);
    protected static final String MATRIX_FLOAT_FORMAT = "%.3f";

    // type of ramp enum - determines data flow direction
    public enum Type {
        OffRamp,
        OnRamp
    }

    // local vars
    protected Type              type;
    protected File              file;
    protected JSONObject        info;

    public RampBase(String filename, Type type) throws IOException  {
        this.type = type;
        this.file = new File(filename);
        logger.info("opening ramp. file: " + file + ", type: " + type);
    }

    protected String getLocFilenameSuffix(Locatable loc) {
        return String.format("%s-%d-%d", loc.getContig(), loc.getStart(), loc.getEnd());
    }

    public void close() throws IOException {
    }

    public Type getType() {
        return type;
    }
}
