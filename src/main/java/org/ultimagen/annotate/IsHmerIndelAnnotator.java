package org.ultimagen.annotate;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.List;
import java.util.stream.Collectors;

public class IsHmerIndelAnnotator extends AnnotatorBase {

    public IsHmerIndelAnnotator(final ReferenceSequenceFile referenceSequenceFile) {
        super(referenceSequenceFile);
    }

    // "hmer_indel_length" and "hmer_indel_nuc"
    public void annotate(final VariantContext vc) {

        String      indelClass = vc.getAttributeAsString(A_INDEL_CLASSIFY, null);
        if ( vc.isIndel() ) {
            List<Allele> alt = vc.getAlleles().stream()
                    .filter(allele -> !allele.isReference())
                    .collect(Collectors.toList());
            if ( alt.size() == 1 ) {
                Allele      altAllele = alt.get(0);
                if (C_INSERT.equals(indelClass)) {
                    /*
                    alt = [x for x in rec['alleles'] if x != rec['ref']][0][1:]
                    if len(set(alt)) != 1:
                        return (0, None)
                    elif fasta_idx[rec['chrom']][rec['pos']].seq.upper() != alt[0]:

                        return (0, None)
                    else:
                        return (utils.hmer_length(fasta_idx[rec['chrom']], rec['pos']), alt[0])
                     */
                    if (getReferenceNucleoid(vc.getContig(), vc.getStart()) == altAllele.getBases()[0]) {
                        annotate(vc, A_HMER_INDEL_LENGTH, hmerLength(vc.getContig(), vc.getStart()));
                        annotate(vc, A_HMER_INDEL_NUC, Character.toString((char) altAllele.getBases()[0]));
                    }

                } else if (C_DELETE.equals(indelClass)) {
                    /*
                    del_seq = rec['ref'][1:]
                    if len(set(del_seq)) != 1:
                        return (0, None)
                    elif fasta_idx[rec['chrom']][rec['pos'] + len(rec['ref']) - 1].seq.upper() != del_seq[0]:
                        return (0, None)
                    else:
                        return (len(del_seq) + utils.hmer_length(fasta_idx[rec['chrom']],
                     */
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
}
