package org.ultimagen.annotate;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.List;
import java.util.stream.Collectors;

public class ConcordanceAnnotator {

    // annotation names
    private final String   A_INDEL_CLASSIFY = "indel_classify";
    private final String   A_INDEL_LENGTH = "indel_length";
    private final String   A_HMER_INDEL_LENGTH = "hmer_indel_length";
    private final String   A_HMER_INDEL_NUC = "hmer_indel_nuc";

    // additional constants
    private final String   C_INSERT = "ins";
    private final String   C_DELETE = "del";

    // bean fields
    private ReferenceSequenceFile referenceSequenceFile;

    public ConcordanceAnnotator(ReferenceSequenceFile referenceSequenceFile) {
        this.referenceSequenceFile = referenceSequenceFile;
    }

    // "indel_classify" and "indel_length"
    public void annotateClassifyIndel(final VariantContext vc) {

        if ( vc.isIndel() ) {

            final int maxAlleleLength = vc.getAlleles().stream()
                    .filter(allele -> !allele.isReference())
                    .map(allele -> allele.length())
                    .max(Integer::compare).get();
            annotate(vc, A_INDEL_CLASSIFY, (vc.getReference().length() < maxAlleleLength) ? C_INSERT : C_DELETE);

            final int refLength = vc.getReference().length();
            annotate(vc, A_INDEL_LENGTH, vc.getAlleles().stream()
                    .filter(allele -> !allele.isReference())
                    .map(allele -> Math.abs(refLength - allele.length()))
                    .max(Integer::compare).get());
        }
    }

    // "hmer_indel_length" and "hmer_indel_nuc"
    public void annotateIsHmerIndel(final VariantContext vc) {

        String      indelClass = vc.getAttributeAsString(A_INDEL_CLASSIFY, null);
        if ( vc.isIndel() ) {
            List<Allele>        alt = vc.getAlleles().stream()
                    .filter(allele -> !allele.isReference())
                    .collect(Collectors.toList());
            if ( alt.size() == 1 ) {
                Allele      altAllele = alt.get(0);
                if (C_INSERT.equals(indelClass)) {
                    if (getReferenceNucleoid(vc.getContig(), vc.getStart()) == altAllele.getBases()[0]) {
                        annotate(vc, A_HMER_INDEL_LENGTH, hmerLength(vc.getContig(), vc.getStart()));
                        annotate(vc, A_HMER_INDEL_NUC, Character.toString((char) altAllele.getBases()[0]));
                    }

                } else if (C_DELETE.equals(indelClass)) {
                    if (getReferenceNucleoid(vc.getContig(), vc.getStart() + vc.getReference().length() - 1) == altAllele.getBases()[0]) {
                        annotate(vc, A_HMER_INDEL_LENGTH, altAllele.length() + hmerLength(vc.getContig(), vc.getStart() + vc.getReference().length() - 1));
                        annotate(vc, A_HMER_INDEL_NUC, Character.toString((char) altAllele.getBases()[0]));
                    }
                }
            }
        }
    }

    private int hmerLength(final String contig, final int start) {
        byte        base = getReferenceNucleoid(contig, start);
        int         length = 1;
        while ( getReferenceNucleoid(contig, start + length) == base )
            length++;
        return length;

    }

    // get nucleoid from reference
    private byte getReferenceNucleoid(final String contig, final int start) {
        return referenceSequenceFile.getSubsequenceAt(contig, start, start).getBases()[0];
    }

    // add a single 'cooked' annotation
    private void annotate(final VariantContext vc, final String name, final Object value) {
        if ( value != null ) {
            vc.getAttributes().put(name, value);
        }
    }
}
