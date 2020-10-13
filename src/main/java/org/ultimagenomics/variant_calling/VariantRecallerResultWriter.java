package org.ultimagenomics.variant_calling;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.Tuple;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.ultimagenomics.flow_based_read.read.FlowBasedRead;
import org.ultimagenomics.haplotype_calling.LHWRefView;
import shaded.cloud_nio.com.google.errorprone.annotations.Var;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class VariantRecallerResultWriter {
    PrintWriter     pw;
    boolean         first = true;
    boolean         debugFormat = false;

    public VariantRecallerResultWriter(File file) {
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
                    if ( vc.getType() != VariantContext.Type.MIXED && (vc.getEnd() != vc.getStart()) )
                        pw.print("-" + vc.getEnd());
                    pw.print(" " + loc);
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
                SimpleInterval      vcSpan = new SimpleInterval(vc.getContig(), vc.getStart(), vc.getEnd());
                List<Tuple<Double,String>>      vcLines = new LinkedList<>();
                for ( int evidenceIndex = 0; evidenceIndex < matrix.evidenceCount() ; evidenceIndex++ ) {

                    // determine matrix values
                    boolean         allValuesNegativeInfinity = true;
                    double          sortKey = Double.NEGATIVE_INFINITY;
                    for (int alleleIndex = 0; alleleIndex < matrix.numberOfAlleles(); alleleIndex++) {
                        lineValues[alleleIndex] = values[alleleIndex][evidenceIndex];
                        if ( lineValues[alleleIndex] != Double.NEGATIVE_INFINITY )
                            allValuesNegativeInfinity = false;
                        sortKey = lineValues[alleleIndex];
                    }

                    // lines which have all values of -Inf are complete alignment failures. Ignore them
                    if ( allValuesNegativeInfinity )
                        continue;

                    // determine length in key space
                    GATKRead        read = matrix.evidence().get(evidenceIndex);
                    int         keyspaceLength = 0;
                    if ( read instanceof FlowBasedRead )
                        keyspaceLength = ((FlowBasedRead)read).getKeyLength();

                    // build basic matrix line
                    String line = String.format("%s %d %d %d %d %s",
                            read.getName(),
                            keyspaceLength,
                            read.isDuplicate() ? 1 : 0,
                            read.isReverseStrand() ? 1 : 0,
                            read.getMappingQuality(),
                            StringUtils.join(ArrayUtils.toObject(lineValues), " "));

                    // add bytes at variant location?
                    StringBuilder bases = new StringBuilder();
                    int           firstBaseUnclippedOfs = 0;
                    if ( read.getContig() != null ) {
                        SimpleInterval  readSpan = new SimpleInterval(read.getContig(), read.getStart(), read.getEnd());
                        if ( readSpan.contains(vcSpan) ) {
                            int ofs = vcSpan.getStart() - readSpan.getStart();
                            int vcLength = vcSpan.getEnd() - vc.getStart() + 1;
                            for (int i = 0; i < vcLength; i++) {
                                int readOfs = getOffsetOnRead(read, ofs + i);
                                if (readOfs >= 0) {
                                    bases.append((char) read.getBase(readOfs));
                                    if ( !read.isReverseStrand() )
                                        firstBaseUnclippedOfs = readOfs + (read.getStart() - read.getUnclippedStart());
                                    else
                                        firstBaseUnclippedOfs = readOfs;
                                } else {
                                    // we don't like '?' bases anymore!
                                    bases.setLength(0);
                                    break;
                                }
                            }
                        }
                    }

                    if ( bases.length() > 0 ) {
                        line += String.format(" %s %d %s", bases, firstBaseUnclippedOfs, read.getReadGroup());
                        vcLines.add(new Tuple<>(sortKey, line));
                    }
                }

                // optional: sort vcLines on second column if present
                vcLines.sort((o1, o2) -> {
                    return -Double.compare(o1.a, o2.a);
                });

                // pour into output file
                vcLines.forEach(doubleStringTuple -> pw.println(doubleStringTuple.b));
            }
        });
    }

    private int getOffsetOnRead(GATKRead read, int ofs) {
        int     readOfs = 0;

        Iterator<CigarElement> iter = read.getCigar().iterator();
        while ( iter.hasNext() ) {
            CigarElement    elem = iter.next();
            CigarOperator   op = elem.getOperator();
            if ( op.consumesReadBases() ) {
                if (ofs < elem.getLength() )
                    return readOfs + ofs;
                else
                    readOfs += elem.getLength();
            }
            ofs -= (op.consumesReferenceBases() ? elem.getLength() : 0);
        }

        // if here, not found
        return -1;
    }
}
