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
    long    umiQualityDropped;

    long    startedAt = System.currentTimeMillis();

    public void writeJson(final File f) {

        final JSONObject  obj = new JSONObject();

        // insert fields
        obj.put("readsIn", readsIn);
        obj.put("bpIn", bpIn);
        putPercentage(obj, "readsOut", readsOut, readsIn);
        putPercentage(obj, "bpOut", bpOut, bpIn);

        putPercentage(obj, "adapter5p", adapter5p, readsIn);
        putPercentage(obj, "adapterMiddle", adapterMiddle, readsIn);
        putPercentage(obj, "adapter3p", adapter3p, readsIn);

        putPercentage(obj, "rqDropped", rqDropped, readsIn);
        putPercentage(obj, "read1TooShortDropped", read1TooShortDropped, readsIn);
        putPercentage(obj, "read2TooShortDropped", read2TooShortDropped, readsIn);
        putPercentage(obj, "bpCutoff", bpCutoff, bpIn);
        putPercentage(obj, "umiQualityDropped", umiQualityDropped, readsIn);

        // performance
        long elapsedMillis = System.currentTimeMillis() - startedAt;
        if ( readsIn != 0 ) {
            obj.put("readsInPerMinute", readsIn / (elapsedMillis / 60000.0f));
        }

        // write
        try (final Writer writer = new PrintWriter(new FileOutputStream(f))) {
            obj.write(writer, 1, 1);
        } catch (IOException e) {
            throw new GATKException("failed to open for writing: " + f, e);
        }
    }

    private void putPercentage(JSONObject obj, String name, long value, long total) {
        StringBuilder       sb = new StringBuilder();
        sb.append(value);
        if ( total != 0 ) {
            sb.append(String.format(" (%.2f%%)", 100.0f * value / total));
        }
        obj.put(name, sb.toString());
    }
}
