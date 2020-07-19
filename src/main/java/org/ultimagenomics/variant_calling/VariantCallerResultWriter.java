package org.ultimagenomics.variant_calling;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.ultimagenomics.haplotype_calling.LHWRefView;
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
    boolean         first = true;
    boolean         debugFormat = false;

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

    public void add(Locatable loc, Map<Integer, AlleleLikelihoods<GATKRead, Allele>> genotypeLikelihoods, List<VariantContext> variants, AssemblyResultSet assemblyResult) {

        // build a map of vcs by startPos
        Map<Integer, VariantContext>    vcStartPos = new LinkedHashMap<>();
        variants.forEach(vc -> {
            vcStartPos.put(vc.getStart(), vc);
        });

        // print location (as a separator)
        if ( debugFormat ) {
            if (first)
                first = false;
            else
                pw.println("");
            pw.println("loc: " + loc);
            pw.println("ref: " + new String(assemblyResult.getFullReferenceWithPadding()));
        }

        // loop on result
        genotypeLikelihoods.forEach((startPos, likelihoods) -> {

            // DK: map to vc? ignore unmapped?
            VariantContext      vc = vcStartPos.get(startPos);
            if ( vc != null ) {

                if ( debugFormat ) {
                    pw.println("");
                    pw.println("variant: " + vc.getContig() + ":" + vc.getStart());
                    pw.println("variant-info: " + vc);

                    // reads
                    pw.println("");
                    pw.println("reads: " + likelihoods.evidenceCount());
                    likelihoods.sampleEvidence(0).forEach(read -> {
                        pw.println("read: " + read);
                    });

                    // alleles
                    pw.println("");
                    pw.println("alleles: " + likelihoods.alleles().size());
                    likelihoods.alleles().forEach(allele -> {
                        pw.println("allele: " + allele);
                    });
                } else {
                    pw.print("#" + vc.getContig() + ":" + vc.getStart());
                    likelihoods.alleles().forEach(allele -> {
                        pw.print(" " + allele);
                    });
                    pw.println("");
                }

                // matrix
                if ( debugFormat ) {
                    pw.println("");
                    pw.println("matrix:");
                }
                LikelihoodMatrix<GATKRead, Allele> matrix = likelihoods.sampleMatrix(0);
                double[][] values = new double[matrix.numberOfAlleles()][matrix.evidenceCount()];
                for (int alleleIndex = 0; alleleIndex < matrix.numberOfAlleles(); alleleIndex++)
                    matrix.copyAlleleLikelihoods(alleleIndex, values[alleleIndex], 0);
                double[] lineValues = new double[matrix.numberOfAlleles()];
                for ( int evidenceIndex = 0; evidenceIndex < matrix.evidenceCount() ; evidenceIndex++ ) {

                    for (int alleleIndex = 0; alleleIndex < matrix.numberOfAlleles(); alleleIndex++)
                        lineValues[alleleIndex] = values[alleleIndex][evidenceIndex];

                    String line = StringUtils.join(ArrayUtils.toObject(lineValues), " ");

                    pw.println(matrix.evidence().get(evidenceIndex).getName() + " " + line);
                }
            }
        });
    }
}
