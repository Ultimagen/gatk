package org.ultimagenomics.ramps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
import java.util.Enumeration;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
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
    private File        file;

    // zip streams
    private JSONObject          info;
    private ZipOutputStream     outputZip;
    private ZipFile             inputZip;

    public OnOffRamp(String filename, Type type) throws IOException  {
        this.type = type;
        this.file = new File(filename);
        logger.info("opening ramp. file: " + file + ", type: " + type);

        // open stream
        if ( type == Type.OffRamp ) {

            // open zip for writing
            this.file.getParentFile().mkdirs();
            outputZip = new ZipOutputStream(new FileOutputStream(this.file));

            // create info object
            info = new JSONObject();
            info.put("created", ZonedDateTime.now(ZoneOffset.UTC).format(DateTimeFormatter.ISO_INSTANT));
            info.put("regions", new JSONArray());

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

    private String getLocFilenameSuffix(Locatable loc) {
        return String.format("%s-%d-%d", loc.getContig(), loc.getStart(), loc.getEnd());
    }

    public void close() throws IOException {

        if ( type == Type.OffRamp ) {

            // add info
            addEntry(null,"info.json", info.toString(2).getBytes());

            // close file
            outputZip.close();

        } else if ( type == Type.OnRamp ) {

            // close file
            inputZip.close();
        }

        logger.info("closing ramp. file: " + file);
    }

    public synchronized  void add(Locatable loc, String name, AlleleLikelihoods<GATKRead, Haplotype> value) throws IOException {

        // add region
        JSONObject regionObj = new JSONObject();
        regionObj.put("contig", loc.getContig());
        regionObj.put("start", loc.getStart());
        regionObj.put("end", loc.getEnd());

        // add global info
        addInfo(regionObj , name + ".haplotypeCount", value.numberOfAlleles());
        addInfo(regionObj, name + ".readCount", value.evidenceCount());

        // add haplotypes
        addHaplotypes(loc,name + ".haplotypes", value.alleles());

        // loop on samples
        for ( int sampleIndex = 0 ; sampleIndex < value.numberOfSamples() ; sampleIndex++ ) {

            // establish context
            String                                 baseName = name + ".samples." + value.getSample(sampleIndex);
            LikelihoodMatrix<GATKRead, Haplotype>  sampleMatrix = value.sampleMatrix(sampleIndex);

            // add to info
            JSONObject      info2 = new JSONObject();
            info2.put("readCount", sampleMatrix.evidenceCount());
            addInfo(regionObj, baseName, info2);

            // write reads
            addReads(loc,baseName + ".reads", sampleMatrix.evidence());

            // write matrix itself
            addSampleMatrix(loc, baseName + ".matrix", sampleMatrix);
        }

        // add to regions
        info.getJSONArray("regions").put(regionObj);
    }

    private void addSampleMatrix(Locatable loc, String name, LikelihoodMatrix<GATKRead, Haplotype> value) throws IOException {

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
        addEntry(loc, name, os.toByteArray());
    }

    private void addHaplotypes(Locatable loc, String name, List<Haplotype> value) throws IOException {

        // build text representation
        ByteArrayOutputStream   os = new ByteArrayOutputStream();
        PrintWriter             pw = new PrintWriter(os);
        pw.println("contig,start,end,ref,cigar,bases");
        for ( Haplotype haplotype : value ) {
            Locatable        hapLoc = haplotype.getGenomeLocation();
            pw.println(String.format("%s,%d,%d,%d,%s,%s",
                    hapLoc.getContig(), hapLoc.getStart(), hapLoc.getEnd(),
                    haplotype.isReference() ? 1 : 0,
                    haplotype.getCigar().toString(),
                    new String(haplotype.getBases())
                    )
            );
        }
        pw.close();

        // add entry
        addEntry(loc, name, os.toByteArray());
    }

    public void addReads(Locatable loc, String name, List<GATKRead> value) throws IOException {

        // build text representation
        ByteArrayOutputStream   os = new ByteArrayOutputStream();
        PrintWriter             pw = new PrintWriter(os);
        pw.println("name,contig,unc_start,start,end,unc_end,cigar");
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
        addEntry(loc, name, os.toByteArray());
    }

    public void add(Locatable loc, String name, Object value) throws IOException {

        // this is a useless default
        addEntry(loc, name, value.toString().getBytes());
    }

    private void addInfo(JSONObject info, String name, Object obj) {

        JSONObject      parent = info;
        String[]        toks =  name.split("\\.");
        for ( int tokIndex = 0 ; tokIndex < toks.length - 1 ; tokIndex++ ) {
            if (!parent.has(toks[tokIndex]))
                parent.put(toks[tokIndex], new JSONObject());
            parent = parent.getJSONObject(toks[tokIndex]);
        }
        parent.put(toks[toks.length - 1], obj);
    }

    private void addEntry(Locatable loc, String name, byte[] bytes) throws IOException {

        String      prefix = loc != null ? getLocFilenameSuffix(loc) + "/" : "";
        ZipEntry    e = new ZipEntry(prefix + name);
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

    public Type getType() {
        return type;
    }



}
