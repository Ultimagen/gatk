package org.ultimagenomics.variant_calling;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import shaded.cloud_nio.com.google.errorprone.annotations.Var;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class VariantCallerResultWriter {
    PrintWriter     pw;

    public VariantCallerResultWriter(File file) {
        try {
            pw = new PrintWriter(file);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public void close() {
        pw.close();
        pw = null;
    }

    public void add(Locatable loc, Map<Integer, AlleleLikelihoods<GATKRead, Allele>> genotypeLikelihoods, List<VariantContext> variants) {

        // build a map of vcs by startPos
        Map<Integer, VariantContext>    vcStartPos = new LinkedHashMap<>();
        variants.forEach(vc -> {
            vcStartPos.put(vc.getStart(), vc);
        });

        // print location (as a separator)
        pw.println("loc: " + loc);

        // loop on result
        genotypeLikelihoods.forEach((startPos, likelihoods) -> {

            // DK: map to vc? ignore unmapped?
            VariantContext      vc = vcStartPos.get(startPos);
            if ( vc != null ) {

                pw.println("variant: " + vc.getContig() + ":" + vc.getStart());
                pw.println("variant-info: " + vc);

                // reads
                pw.println("reads: " + likelihoods.evidenceCount());
                likelihoods.sampleEvidence(0).forEach(read -> {
                    pw.println("read: " + read);
                });

                // alleles
                pw.println("alleles: " + likelihoods.alleles().size());
                likelihoods.alleles().forEach(allele -> {
                    pw.println("allele: " + allele);
                });

                // matrix
                pw.println("matrix:");
                LikelihoodMatrix<GATKRead, Allele> matrix = likelihoods.sampleMatrix(0);
                double[] values = new double[matrix.evidenceCount()];
                for (int alleleIndex = 0; alleleIndex < matrix.numberOfAlleles(); alleleIndex++) {

                    matrix.copyAlleleLikelihoods(alleleIndex, values, 0);
                    String line = StringUtils.join(ArrayUtils.toObject(values), ",");

                    pw.println(line);
                }
            }
        });
    }
}
