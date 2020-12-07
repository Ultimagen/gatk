package org.ultimagenomics.ramps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import htsjdk.samtools.util.Locatable;
import org.json.JSONObject;
import org.json.JSONTokener;

import java.io.*;
import java.time.ZoneOffset;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Enumeration;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

public class OnOffRamp {

    private static final Logger logger = LogManager.getLogger(OnOffRamp.class);

    // type of ramp enum - determines data flow direction
    public enum Type {
        OffRamp,
        OnRamp
    }

    // local vars
    private Type        type;
    private Locatable   loc;
    private File        file;

    // zip streams
    private JSONObject          info;
    private ZipOutputStream     outputZip;
    private ZipFile             inputZip;

    public OnOffRamp(String filename, Locatable loc, Type type) throws IOException  {
        this.type = type;
        this.loc = loc;
        this.file = new File(filename + getLocFilenameSuffix());
        logger.info("opening ramp. file: " + file + ", type: " + type);

        // open stream
        if ( type == Type.OffRamp ) {

            // open zip for writing
            this.file.getParentFile().mkdirs();
            outputZip = new ZipOutputStream(new FileOutputStream(this.file));

            // create info object
            info = new JSONObject();
            info.put("created", ZonedDateTime.now(ZoneOffset.UTC).format(DateTimeFormatter.ISO_INSTANT));
            JSONObject regionObj = new JSONObject();
            regionObj.put("contig", loc.getContig());
            regionObj.put("start", loc.getStart());
            regionObj.put("end", loc.getEnd());
            info.put("region", regionObj);

        } else if ( type == Type.OnRamp ) {

            // open zip for reading
            inputZip = new ZipFile(this.file);

            // read info object
            InputStream     is = getEntry("info.json");
            info = new JSONObject(new JSONTokener(is));

        } else {
            throw new Error("type not supprted: " + type);
        }
    }

    private String getLocFilenameSuffix() {
        return String.format("-%s-%d-%d.zip", loc.getContig(), loc.getStart(), loc.getEnd());
    }

    public void close() throws IOException {

        if ( type == Type.OffRamp ) {

            // add info
            addEntry("info.json", info.toString(2).getBytes());

            // close file
            outputZip.close();

        } else if ( type == Type.OnRamp ) {

            // close file
            inputZip.close();
        }

        logger.info("closing ramp. file: " + file);
    }

    public void add(String name, AlleleLikelihoods<GATKRead, Haplotype> value) throws IOException {

        // add global info
        JSONObject      info = new JSONObject();
        info.put("haplotypeCount", value.numberOfAlleles());
        info.put("readCount", value.evidenceCount());
        addInfo(name, info);

        // add haplotypes
        addHaplotypes(name + ".haplotypes", value.alleles());

        // loop on samples
        for ( int sampleIndex = 0 ; sampleIndex < value.numberOfSamples() ; sampleIndex++ ) {

            // establish context
            String                                 baseName = name + ".samples." + value.getSample(sampleIndex);
            LikelihoodMatrix<GATKRead, Haplotype>  sampleMatrix = value.sampleMatrix(sampleIndex);

            // add to info
            JSONObject      info1 = new JSONObject();
            info1.put("readCount", sampleMatrix.evidenceCount());
            addInfo(baseName, info1);

            // write reads
            addReads(baseName + ".reads", sampleMatrix.evidence());

            // write matrix itself
            addSampleMatrix(baseName + ".matrix", sampleMatrix);
        }
    }

    private void addSampleMatrix(String name, LikelihoodMatrix<GATKRead, Haplotype> value) throws IOException {

        // build text representation
        ByteArrayOutputStream   os = new ByteArrayOutputStream();
        PrintWriter             pw = new PrintWriter(os);

        // header line
        pw.print("read");
        for ( int n = 0 ; n < value.numberOfAlleles() ; n++ )
            pw.print(",h" + n);
        pw.println("");

        // walk the lines
        for ( int readIndex  = 0 ; readIndex < value.evidenceCount() ; readIndex++ ) {
            pw.print(value.getEvidence(readIndex).getName());
            for ( int hapIndex = 0 ; hapIndex < value.numberOfAlleles() ; hapIndex++ )
                pw.print("," + value.get(hapIndex, readIndex));
            pw.println("");
        }
        pw.close();

        // add entry
        addEntry(name, os.toByteArray());
    }

    private void addHaplotypes(String name, List<Haplotype> value) throws IOException {

        // build text representation
        ByteArrayOutputStream   os = new ByteArrayOutputStream();
        PrintWriter             pw = new PrintWriter(os);
        pw.println("contig,start,end,ref,cigar,bases");
        for ( Haplotype haplotype : value ) {
            Locatable        loc = haplotype.getGenomeLocation();
            pw.println(String.format("%s,%d,%d,%d,%s,%s",
                    loc.getContig(), loc.getStart(), loc.getEnd(),
                    haplotype.isReference() ? 1 : 0,
                    haplotype.getCigar().toString(),
                    new String(haplotype.getBases())
                    )
            );
        }
        pw.close();

        // add entry
        addEntry(name, os.toByteArray());
    }

    public void addReads(String name, List<GATKRead> value) throws IOException {

        // build text representation
        ByteArrayOutputStream   os = new ByteArrayOutputStream();
        PrintWriter             pw = new PrintWriter(os);
        pw.println("name,contig,start,end");
        for ( GATKRead read : value ) {
            pw.println(String.format("%s,%s,%d,%d,%d,%d,%s",
                    read.getName(),
                    read.getContig(),
                    read.getUnclippedEnd(), read.getStart(),
                    read.getEnd(), read.getUnclippedEnd(),
                    read.getCigar().toString()
                    )
            );
        }
        pw.close();

        // add entry
        addEntry(name, os.toByteArray());
    }

    public void add(String name, Object value) throws IOException {

        // this is a useless default
        addEntry(name, value.toString().getBytes());
    }

    private void addInfo(String name, JSONObject obj) {

        JSONObject      parent = info;
        String[]        toks =  name.split("\\.");
        for ( int tokIndex = 0 ; tokIndex < toks.length - 1 ; tokIndex++ ) {
            if (!parent.has(toks[tokIndex]))
                parent.put(toks[tokIndex], new JSONObject());
            parent = parent.getJSONObject(toks[tokIndex]);
        }
        parent.put(toks[toks.length - 1], obj);
    }

    private void addEntry(String name, byte[] bytes) throws IOException {

        ZipEntry    e = new ZipEntry(name);
        outputZip.putNextEntry(e);
        outputZip.write(bytes);
        outputZip.closeEntry();
    }

    private InputStream getEntry(String name) throws IOException {

        Enumeration<? extends ZipEntry> entries = inputZip.entries();
        while ( entries.hasMoreElements() ) {
            ZipEntry entry = entries.nextElement();
            if ( entry.getName().equals(name) )
                return inputZip.getInputStream(entry);
        }

        // if here, not found
        throw new IOException("no such: " + name);
    }

    public AlleleLikelihoods<GATKRead, Haplotype> getAlleleLikelihoods(String readLikelihoods) {
        return null;
    }

}
