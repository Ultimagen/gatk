package org.ultimagenomics.flow_based_read.tests;

import org.ultimagenomics.flow_based_read.read.FlowBasedRead;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;

public class test_flow_based_read {
    public static void main(String[] args) throws FileNotFoundException, IOException{
        System.out.println("Input file: " + args[0]);
        System.out.println("Output prefix: " + args[1]);
        System.out.println("Flow order: " + args[2]);
        final SamReader reader = SamReaderFactory.makeDefault().open(new File(args[0]));
        final String flowOrder = args[2];
        int count = 0 ;
        Iterator<SAMRecord> i;
        FileWriter fos;
        FlowBasedRead fbr = null;
        for (i = reader.iterator(), count=0; (i.hasNext()) && (count < 100); count++  ) {
            fbr = new FlowBasedRead(i.next(), flowOrder, 13);
            System.out.println(String.format("> Key length: %d, sequence length %d", fbr.totalKeyBases(), fbr.seqLength()));

            fbr.apply_alignment();
            System.out.println(String.format("< Key length: %d, sequence length %d", fbr.totalKeyBases(), fbr.seqLength()));
            if (fbr.totalKeyBases()!=fbr.seqLength()) {
                System.out.println("!@#@!$#%#%#%$@$@#!");
            }
            fos = new FileWriter(args[1] + "." + Integer.toString(count) + ".key.txt");
            fbr.writeKey(fos);
            fos.close();
            fos = new FileWriter(args[1] + "." + Integer.toString(count) + ".matrix.txt");
            fbr.writeMatrix(fos);
            fos.close();

        }
    }
}
