package org.ultimagenomics.flow_based_read.tests;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.ultimagenomics.flow_based_read.alignment.FlowBasedAlignmentEngine;
import org.ultimagenomics.flow_based_read.read.FlowBasedHaplotype;
import org.ultimagenomics.flow_based_read.read.FlowBasedRead;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.genomicsdb.importer.model.ChromosomeInterval;
import org.ultimagenomics.flow_based_read.utils.FlowBasedAlignmentArgumentCollection;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;


public class test_haplotype_matching {
    public static void main(String[] args) throws IOException{
        System.out.println("Input reads file: " + args[0]);
        System.out.println("Input haplotypes file: " + args[1]);
        System.out.println("CR: " + args[2]);
        System.out.println("Output prefix: " + args[1]);
        final SamReader reader = SamReaderFactory.makeDefault().open(new File(args[0]));
        ArrayList<GATKRead> reads = new ArrayList<>();
        ArrayList<Haplotype> haplotypes = new ArrayList<>();
        int count_reads = 0;
        int count_haplotypes = 0;
        ArrayList<SAMFileHeader> headers = new ArrayList<>();
        for ( final SAMRecord rec : reader ) {
            SAMRecordToGATKReadAdapter tmp = new SAMRecordToGATKReadAdapter(rec);
            if (!tmp.getAttributeAsString("RG").startsWith("ArtificialHaplotype")) {
                if (tmp.getAttributeAsString(("CR")).equals(args[2])) {
                    reads.add(tmp);
                    count_reads++;
                }
            } else {
                if (tmp.getAttributeAsString(("CR")).equals(args[2])) {
                    ChromosomeInterval gl = new ChromosomeInterval(tmp.getContig(), tmp.getStart(), tmp.getEnd());

                    Haplotype hap = new Haplotype(tmp.getBases(), gl);
                    hap.setCigar(tmp.getCigar());
                    haplotypes.add(hap);
                    count_haplotypes++;
                }
            }

        }

        System.out.println(String.format("%d reads %d haplotypes", count_reads, count_haplotypes));
        ArrayList<FlowBasedRead> fbrs = new ArrayList<>();
        int cnt = 0 ;
        FlowBasedAlignmentArgumentCollection fbargs = new FlowBasedAlignmentArgumentCollection();
        for (final GATKRead r :  reads) {
            fbrs.add(new FlowBasedRead(r, "TACG", 8, fbargs));
        }
        for (FlowBasedRead fbr : fbrs ) fbr.apply_alignment();

        ArrayList<FlowBasedHaplotype> fbhs = new ArrayList<>();
        for (final Haplotype hap : haplotypes) {
            fbhs.add(new FlowBasedHaplotype(hap, "TACG", 8));
        }

        FlowBasedAlignmentEngine fbe = new FlowBasedAlignmentEngine(new FlowBasedAlignmentArgumentCollection(), -5, 0.02);

        final AlleleLikelihoods<GATKRead, Haplotype> haplotypeReadLikelihoods =
                fbe.computeReadLikelihoods(haplotypes, reads, true);

    }
}
