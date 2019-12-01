package org.ultimagenomics.flow_based_read.tests;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class AlleleLikelihoodWriter<EVIDENCE extends Locatable, A extends Allele> implements AutoCloseable {
    Path outputPath;
    SimpleInterval outputInterval;
    FileWriter output;
    public AlleleLikelihoodWriter(final Path _outputPath, final SimpleInterval _interval) {
        this.outputPath = _outputPath;
        this.outputInterval = _interval;
        try {
            output = new FileWriter(outputPath.toString());
        } catch (IOException err) {
            throw new IllegalArgumentException(String.format("Unable to open %s for writing", _outputPath.toString()));
        }
    }

    public void writeAlleleLikelihoods(AlleleLikelihoods<GATKRead, Haplotype> likelihoods){
        List<String> samples = likelihoods.samples();
        List<Haplotype> haplotypes = likelihoods.alleles();

        Haplotype first_hap = haplotypes.get(0);
        try {
            if (first_hap.getLocation().getContig().equals(outputInterval.getContig())) {
                if (((first_hap.getStartPosition() >= outputInterval.getStart()) && (first_hap.getStartPosition() <= outputInterval.getEnd())) &&
                        ((first_hap.getStopPosition() >= outputInterval.getStart()) && (first_hap.getStopPosition() <= outputInterval.getEnd()))) {

                    output.write(String.format("> Location %s:%d-%d\n", first_hap.getLocation().getContig(),
                            first_hap.getStartPosition(), first_hap.getStopPosition()));
                    output.write(">> Haplotypes\n");
                    for (int i = 0; i < haplotypes.size(); i++) {
                        output.write(String.format("%04d\t%s\n", i, haplotypes.get(i).toString()));
                    }
                    for (int s = 0; s < samples.size(); s++) {
                        output.write(String.format(">> Sample %s\n", samples.get(s)));
                        output.write(">>> Reads\n");
                        List<GATKRead> reads = likelihoods.sampleEvidence(s);
                        for (int i = 0; i < reads.size(); i++) {
                            output.write(String.format("%04d\t%s\n", i, reads.get(i).getName()));
                        }
                        output.write(">>> Matrix\n");
                        for (int allele = 0; allele < likelihoods.sampleMatrix(s).numberOfAlleles(); allele++) {
                            for (int read = 0; read < likelihoods.sampleMatrix(s).evidenceCount(); read++) {
                                output.write(Double.toString(likelihoods.sampleMatrix(s).get(allele, read)));
                                output.write(" ");
                            }
                            output.write("\n");
                        }
                    }
                }
            }
        } catch (IOException err) {
            throw new RuntimeException(String.format("Unable to write matrix to file"));
        }
    }



    @Override
    public void close() {
        try {
            output.close();
        } catch (IOException err) {
            throw (new RuntimeException("Unable to close matrix file"));
        }
    }
}
