package org.ultimagen.flowBasedRead.tests;
import htsjdk.samtools.*;
import org.ultimagen.flowBasedRead.utils.Variant;
import java.io.File;

/**
 * @author ilya
 * Test class - tests parsing SAM file with uncertainty tags
 */

public class test_sam_parsing {
    public static void main(String args[]) {
        System.out.println("Input file: " + args[0]);
        final SamReader reader = SamReaderFactory.makeDefault().open(new File(args[0]));
        int i = 0 ;
        for (final SAMRecord samRecord : reader ) {
            float likelihood = (float) samRecord.getAttribute("LL");

            String alt = (String) samRecord.getAttribute("AL");
            if (alt.length() == 0) {
                continue;
            }
            else {
                String[] mutations = alt.split(";");
                Variant[] variants = new Variant[mutations.length];
                for (int j = 0; j<mutations.length; j++) {
                    variants[j] = new Variant(mutations[j]);
                    System.out.println(String.format("%.2f %s",
                            likelihood,
                            variants[j].toString()));
                }
            }
        }
        System.out.println(Integer.toString(i) + " reads");
    }

    //private static
}

