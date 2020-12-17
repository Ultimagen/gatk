package org.ultimagenomics.ramps;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import htsjdk.samtools.util.Locatable;
import org.json.JSONArray;
import org.json.JSONObject;
import org.json.JSONTokener;

import java.io.*;
import java.time.ZoneOffset;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

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
