package org.ultimagenomics.flow_based_read.tests;

import org.ultimagenomics.flow_based_read.read.FlowBasedRead;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.ultimagenomics.flow_based_read.utils.FlowBasedAlignmentArgumentCollection;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;

public class test_flow_based_read {
    public static void main(String[] args) throws FileNotFoundException, IOException{
        FlowBasedAlignmentArgumentCollection fbargs = new FlowBasedAlignmentArgumentCollection();

        System.out.println("Input file: " + args[0]);
        System.out.println("Output prefix: " + args[1]);
        System.out.println("Flow order: " + args[2]);

        if (args.length > 4 && args[4].equals("simulate")) {
            fbargs.probability_ratio_threshold = 0.01;
//            fbargs.remove_longer_than_one_indels = true;
            fbargs.lump_probs = true;
            fbargs.only_ins_or_del = true;
            fbargs.remove_one_to_zero_probs = true;
        }
        final SamReader reader = SamReaderFactory.makeDefault().open(new File(args[0]));
        final String flowOrder = args[2];
        final int limitCount = (args.length > 3) ? Integer.valueOf(args[3]) : 100;
        System.out.println("limitCount: " + limitCount);
        int count = 0 ;
        Iterator<SAMRecord> i;
        FileWriter fos;
        FlowBasedRead fbr = null;
        for (i = reader.iterator(), count=0; (i.hasNext()) && (count < limitCount); count++  ) {
            fbr = new FlowBasedRead(i.next(), flowOrder, 12, fbargs);
            if ( limitCount < 1000 )
                System.out.println(String.format("> Key length: %d, sequence length %d", fbr.totalKeyBases(), fbr.seqLength()));

            fbr.apply_alignment();
            if ( limitCount < 1000 )
                System.out.println(String.format("< Key length: %d, sequence length %d", fbr.totalKeyBases(), fbr.seqLength()));
            if (fbr.totalKeyBases()!=fbr.seqLength()) {
                System.out.println("!@#@!$#%#%#%$@$@#!");
            }
            if ( limitCount < 1000 ) {
                fos = new FileWriter(args[1] + "." + Integer.toString(count) + ".key.txt");
                fbr.writeKey(fos);
                fos.close();
                fos = new FileWriter(args[1] + "." + Integer.toString(count) + ".matrix.txt");
                fbr.writeMatrix(fos);
                fos.close();
            } else if ( (count % 1000) == 0 ) {
                System.out.println("record count: " + count);

            }

        }
        System.out.println("total records read: " + count);
    }
}
