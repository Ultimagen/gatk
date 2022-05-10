package org.broadinstitute.hellbender.tools.walkers.pipeline;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.json.JSONObject;

import java.io.*;

public class SingleCellPipelineToolStatistics {

    long    readsIn;
    long    readsOut;
    long    bpIn;
    long    bpOut;

    long    adapter5p;
    long    adapterMiddle;
    long    adapter3p;

    long    rqDropped;
    long    read1TooShortDropped;
    long    read2TooShortDropped;
    long    bpCutoff;

    public void writeJson(final File f) {

        final JSONObject  obj = new JSONObject();

        obj.put("readsIn", readsIn);
        obj.put("readsOut", readsOut);
        obj.put("bpIn", bpIn);
        obj.put("bpOut", bpOut);

        obj.put("adapter5p", adapter5p);
        obj.put("adapterMiddle", adapterMiddle);
        obj.put("adapter3p", adapter3p);

        obj.put("rqDropped", rqDropped);
        obj.put("read1TooShortDropped", read1TooShortDropped);
        obj.put("read2TooShortDropped", read2TooShortDropped);
        obj.put("bpCutoff", bpCutoff);

        try (final Writer writer = new PrintWriter(new FileOutputStream(f))) {
            obj.write(writer, 1, 1);
        } catch (IOException e) {
            throw new GATKException("failed to open for writing: " + f, e);
        }
    }
}
