package org.ultimagenomics.ramps;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.ByteArrayOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.ZoneOffset;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

public class PreFilterOffRamp extends RampBase {

    private ZipOutputStream outputZip;

    public PreFilterOffRamp(String filename) throws IOException {
        super(filename, Type.OffRamp);

        // open zip for writing
        this.file.getParentFile().mkdirs();
        outputZip = new ZipOutputStream(new FileOutputStream(this.file));

        // create info object
        info = new JSONObject();
        info.put("created", ZonedDateTime.now(ZoneOffset.UTC).format(DateTimeFormatter.ISO_INSTANT));
        info.put("regions", new JSONArray());
    }

    public void close() throws IOException {

        // add info
        addEntry(null,"info.json", info.toString(2).getBytes());

        // close file
        outputZip.close();

        super.close();
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
            LikelihoodMatrix<GATKRead, Haplotype> sampleMatrix = value.sampleMatrix(sampleIndex);

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
        ByteArrayOutputStream os = new ByteArrayOutputStream();
        PrintWriter pw = new PrintWriter(os);

        // header line
        pw.print("read");
        for ( int n = 0 ; n < value.numberOfAlleles() ; n++ )
            pw.print(",h" + n);
        pw.println("");

        // walk the lines
        for ( int readIndex  = 0 ; readIndex < value.evidenceCount() ; readIndex++ ) {
            pw.print(value.getEvidence(readIndex).getName());
            for ( int hapIndex = 0 ; hapIndex < value.numberOfAlleles() ; hapIndex++ )
                pw.print("," + String.format(MATRIX_FLOAT_FORMAT, value.get(hapIndex, readIndex)));
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
        pw.println("contig,start,end,ref,cigar,bases,score,alignmentStartHapwrtRef");
        for ( Haplotype haplotype : value ) {
            Locatable        hapLoc = haplotype.getGenomeLocation();
            pw.println(String.format("%s,%d,%d,%d,%s,%s,%.5f,%d",
                    hapLoc.getContig(), hapLoc.getStart(), hapLoc.getEnd(),
                    haplotype.isReference() ? 1 : 0,
                    haplotype.getCigar().toString(),
                    new String(haplotype.getBases()),
                    haplotype.getScore(),
                    haplotype.getAlignmentStartHapwrtRef()
                    )
            );
        }
        pw.close();

        // add entry
        addEntry(loc, name, os.toByteArray());
    }

    public void addReads(Locatable loc, String name, List<GATKRead> value) throws IOException {

        /**
         * at this point it is not clear if we need the reads at all. we'll keep a skeleton here
         */

        // build text representation
        ByteArrayOutputStream   os = new ByteArrayOutputStream();
        PrintWriter             pw = new PrintWriter(os);
        pw.println("name");
        for ( GATKRead read : value ) {
            pw.println(String.format("%s",
                    read.getName()
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
        ZipEntry e = new ZipEntry(prefix + name);
        outputZip.putNextEntry(e);
        outputZip.write(bytes);
        outputZip.closeEntry();
    }

}
