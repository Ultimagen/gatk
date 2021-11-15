package org.ultimagen.flowBasedRead.utils;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

/*
A class for logging likelihood matrics, possibly constrained by an (include-only) output interval

The output is textual. Each matrix contains several sections. As section begins with a heading line (staring with >> <heading-name>(
and follows with the section data organised as a two dimentional matrix over several lines.
 */
public class AlleleLikelihoodWriter<EVIDENCE extends Locatable, A extends Allele> implements AutoCloseable {
    final Path outputPath;
    final SimpleInterval outputInterval;
    final FileWriter output;


    public AlleleLikelihoodWriter(final Path _outputPath, final SimpleInterval _interval) {
        this.outputPath = _outputPath;
        this.outputInterval = _interval;
        try {
            output = new FileWriter(outputPath.toString());
        } catch (IOException err) {
            throw new IllegalArgumentException(String.format("Unable to open %s for writing", _outputPath.toString()));
        }
    }

    /**
     * Add a likelihood matrix to the output. Only haplotypes falling within the output interval will be output
     * @param likelihoods - matrix to add
     */
    public void writeAlleleLikelihoods(final AlleleLikelihoods<GATKRead, Haplotype> likelihoods){
        final List<String> samples = likelihoods.samples();
        final List<Haplotype> haplotypes = likelihoods.alleles();

        final Haplotype first_hap = haplotypes.get(0);
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
                        final List<GATKRead> reads = likelihoods.sampleEvidence(s);
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

                        output.write(">>> Read->Haplotype in Full\n");
                        for (int allele = 0; allele < likelihoods.sampleMatrix(s).numberOfAlleles(); allele++) {
                            for (int read = 0; read < likelihoods.sampleMatrix(s).evidenceCount(); read++) {
                                output.write(String.format("%04d\t%s\t%s\t%04d\t%s\n",
                                        read, reads.get(read).getName(),
                                        Double.toString(likelihoods.sampleMatrix(s).get(allele, read)),
                                        allele, haplotypes.get(allele).toString()));
                            }
                            output.write("\n");
                        }
                    }
                    output.flush();
                }
            }
        } catch (IOException err) {
            throw new RuntimeException(String.format("Unable to write matrix to file"));
        }
    }



    @Override
    public void close() {
        try {
            output.flush();
            output.close();
        } catch (IOException err) {
            throw (new RuntimeException("Unable to close matrix file"));
        }
    }
}
