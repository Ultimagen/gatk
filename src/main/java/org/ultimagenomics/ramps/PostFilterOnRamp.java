package org.ultimagenomics.ramps;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.json.JSONObject;
import org.json.JSONTokener;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

public class PostFilterOnRamp extends RampBase {

    private ZipFile inputZip;

    public PostFilterOnRamp(String filename) throws IOException {
        super(filename, Type.OnRamp);

        // open zip for reading
        inputZip = new ZipFile(this.file);

        // read info object
        InputStream is = getEntry(null,"info.json");
        info = new JSONObject(new JSONTokener(is));
    }

    public void close() throws IOException {

        // close file
        inputZip.close();

        super.close();
    }

    private InputStream getEntry(Locatable loc, String name) throws IOException {

        String      prefix = loc != null ? getLocFilenameSuffix(loc) + "/" : "";
        name = prefix + name;

        Enumeration<? extends ZipEntry> entries = inputZip.entries();
        while ( entries.hasMoreElements() ) {
            ZipEntry entry = entries.nextElement();
            if ( entry.getName().equals(name) )
                return inputZip.getInputStream(entry);
        }

        // if here, not found
        throw new IOException("no such: " + name);
    }

    public AlleleLikelihoods<GATKRead, Haplotype> getAlleleLikelihoods(Locatable loc, String name, SampleList sampleList, Map<String,List<GATKRead>> readsBySample) throws IOException {

        // get haplotypes
        List<Haplotype> haplotypes = getHaplotypes(loc, name + ".haplotypes");

        // collect reads by name
        Map<String, GATKRead> readsByName = new LinkedHashMap<>();
        for ( List<GATKRead> readList : readsBySample.values() )
            for ( GATKRead read : readList )
                readsByName.put(read.getName(), read);

        // walk the samples
        List<List<GATKRead>>        evidenceBySampleList = new LinkedList<>();
        double[][][]                values = new double[sampleList.numberOfSamples()][][];
        for ( int sampleIndex = 0 ; sampleIndex < sampleList.numberOfSamples() ; sampleIndex++ ) {

            // establish context
            String                  sampleName = sampleList.getSample(sampleIndex);
            String                  baseName = name + ".samples." + sampleName;

            // enrich reads
            List<GATKRead>          enrichedReads = enrichReads(loc,baseName + ".reads", readsByName);
            evidenceBySampleList.add(enrichedReads);

            // read matrix
            double[][]              matrix = getSampleMatrix(loc, baseName + ".matrix", haplotypes, enrichedReads);
            values[sampleIndex] = matrix;
        }

        // create returned object
        IndexedAlleleList<Haplotype> allelList = new IndexedAlleleList<>(haplotypes);
        AlleleLikelihoods<GATKRead, Haplotype>      alleleLikelihoods = new AlleleLikelihoods<>(
                                                        allelList, sampleList, evidenceBySampleList, values);
        for ( int sampleIndex = 0 ; sampleIndex < sampleList.numberOfSamples() ; sampleIndex++ )
            alleleLikelihoods.sampleMatrix(sampleIndex);

            return alleleLikelihoods;
    }

    private List<Haplotype> getHaplotypes(Locatable loc, String name) throws IOException {

        // open csv file
        CSVReader           reader = new CSVReader(getEntry(loc, name));
        int                 contigColumn = reader.getColumn("contig");
        int                 startColumn = reader.getColumn("start");
        int                 endColumn = reader.getColumn("end");
        int                 refColumn = reader.getColumn("ref");
        int                 cigarColumn = reader.getColumn("cigar");
        int                 basesColumn = reader.getColumn("bases");
        int                 scoreColumn = reader.getColumn("score");
        int                 alignmentStartHapwrtRefColumn = reader.getColumn("alignmentStartHapwrtRef");

        // read haplotypes
        List<Haplotype>     haplotypes = new LinkedList<>();
        while ( reader.next() ) {

            Haplotype       haplotype = new Haplotype(reader.get(basesColumn).getBytes(), reader.getBoolean(refColumn));
            haplotype.setGenomeLocation(reader.getLocatable(contigColumn, startColumn, endColumn));
            haplotype.setScore(reader.getDouble(scoreColumn));
            haplotype.setCigar(reader.getCigar(cigarColumn));
            haplotype.setAlignmentStartHapwrtRef(reader.getInteger(alignmentStartHapwrtRefColumn));

            /**
             * what to do about contigs? can they be regenerted?
             */
            // haplotype.contigs = h.contigs;

            haplotypes.add(haplotype);

        }

        // close up
        reader.close();

        return haplotypes;
    }

    private List<GATKRead> enrichReads(Locatable loc, String name, Map<String, GATKRead> readsByName) throws IOException {

        // open csv file
        CSVReader           reader = new CSVReader(getEntry(loc, name));
        int                 nameColumn = reader.getColumn("name");

        // read rows, enrich
        List<GATKRead>      enrichedReads = new LinkedList<>();
        while ( reader.next() ) {

            // locate read
            String          readName = reader.get(nameColumn);
            GATKRead        read = readsByName.get(readName);
            if ( read == null )
                throw new IOException("no such read: " + readName);

            // enrich here ...

            // add to list
            enrichedReads.add(read);
        }

        // close up
        reader.close();

        return enrichedReads;
    }

    private double[][] getSampleMatrix(Locatable loc, String name, List<Haplotype> haplotypes, List<GATKRead> reads) throws IOException {

        // open csv file
        CSVReader           reader = new CSVReader(getEntry(loc, name));
        if ( reader.getColumnCount() != (haplotypes.size() + 1) )
            throw new IOException(String.format("column count %d does not match number of haplotypes %d",
                                                                reader.getColumnCount(), haplotypes.size()));
        // allocate matrix
        // for some (unknown) reason the matrix has one more row at the end with value dup to the one before it???
        double[][]          matrix = new double[haplotypes.size()][];
        for ( int i = 0 ; i < haplotypes.size() ; i++ )
            matrix[i] = new double[reads.size() + 1];
        int                 row = 0;
        while ( reader.next() ) {

            // establish read
            String      readName = reader.get(0);
            if ( !readName.equals(reads.get(row).getName()) )
                throw new IOException(String.format("Matrix read %s does not match expected %s",
                                                        readName, reads.get(row).getName()));
            // read data
            for ( int i = 0 ; i < haplotypes.size() ; i++ )
                matrix[i][row] = reader.getDouble(i + 1);

            // advance
            row++;
        }
        if ( row != reads.size() )
            throw new IOException(String.format("%d rows found, but there are %d reads", row, reads.size()));

        // add the extra row
        for ( int i = 0 ; i < haplotypes.size() ; i++ )
            matrix[i][row] = matrix[i][row - 1];

        // return
        return matrix;
    }

    private class CSVReader {

        private BufferedReader reader;
        Map<String,Integer>             columns = new LinkedHashMap<>();
        String[]                        row;

        CSVReader(InputStream is) throws  IOException {

            // open reader
            reader = new BufferedReader(new InputStreamReader(is));

            // read column headings
            String      toks[] = reader.readLine().split(",");
            for ( int n = 0 ; n < toks.length ; n++ )
                columns.put(toks[n], n);
        }

        void close() throws IOException {
            reader.close();
        }

        boolean next() throws IOException {

            String  line = reader.readLine();
            if ( line == null )
                return false;

            row = line.split(",");
            return true;
        }

        int getColumnCount() {
            return columns.size();
        }

        int getColumn(String name) throws IOException {

            Integer         column = columns.get(name);
            if ( column == null )
                throw new IOException("no such column: " + name);

            return column;
        }

        String get(int column) throws IOException {
            return row[column];
        }

        int getInteger(int column) throws IOException {
            return Integer.parseInt(get(column));
        }

        double getDouble(int column) throws IOException {
            return Double.parseDouble(get(column));
        }

        boolean getBoolean(int column) throws IOException {
            return getInteger(column) != 0;
        }

        Cigar getCigar(int column) throws IOException {
            return TextCigarCodec.decode(get(column));
        }

        Locatable getLocatable(int contigColumn, int startColumn, int endColumn) throws IOException {
            return new SimpleInterval(get(contigColumn), getInteger(startColumn), getInteger(endColumn));
        }
    }

}
