package org.ultimagen.annotate;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public class GetMotifAroundAnnotation extends AnnotatorBase {

    public GetMotifAroundAnnotation(final ReferenceSequenceFile referenceSequenceFile) {
        super(referenceSequenceFile);
    }

    private int         size = 5;

    @Override
    public void annotate(final VariantContext vc) {

        String      indelClass = vc.getAttributeAsString(A_INDEL_CLASSIFY, null);
        int         hmerIndelLength = vc.getAttributeAsInt(A_HMER_INDEL_LENGTH, 0);

        if ( indelClass != null && hmerIndelLength > 0 ) {
            annotateMotifAroundHherIndel(vc, hmerIndelLength);
        } else if ( indelClass != null && hmerIndelLength == 0 ) {
            annotateMotifAroundNonHherIndel(vc);
        } else {
            annotateMotifAroundSnp(vc);
        }
    }

    private void annotateMotifAroundHherIndel(final VariantContext vc, final int hmerLength) {
        /*
        chrom = faidx[rec['chrom']]
        pos = rec['pos']
        hmer_length = rec['hmer_indel_length']
        return chrom[pos - size:pos].seq.upper(), chrom[pos + hmer_length:pos + hmer_length + size].seq.upper()
         */
        vc.getAttributes().put(A_LEFT_MOTIF, getReferenceMotif(vc.getContig(), vc.getStart() - size, vc.getStart()));
        vc.getAttributes().put(A_RIGHT_MOTIF, getReferenceMotif(vc.getContig(), vc.getStart() + hmerLength, vc.getStart() + hmerLength + size));
    }

    private void annotateMotifAroundNonHherIndel(final VariantContext vc) {
        /*
        return chrom[pos - size:pos].seq.upper(),\
            chrom[pos + len(rec['ref']) - 1:pos +
                  len(rec['ref']) - 1 + size].seq.upper()
         */
        int         refLength = vc.getReference().length();
        vc.getAttributes().put(A_LEFT_MOTIF, getReferenceMotif(vc.getContig(), vc.getStart() - size, vc.getStart()));
        vc.getAttributes().put(A_RIGHT_MOTIF, getReferenceMotif(vc.getContig(), vc.getStart() + refLength - 1, vc.getStart() + refLength - 1 + size));
    }

    private void annotateMotifAroundSnp(final VariantContext vc) {
        /*
        return chrom[pos - size - 1:pos - 1].seq.upper(), chrom[pos:pos + size].seq.upper()
         */
        vc.getAttributes().put(A_LEFT_MOTIF, getReferenceMotif(vc.getContig(), vc.getStart() - size - 1, vc.getStart() - 1));
        vc.getAttributes().put(A_RIGHT_MOTIF, getReferenceMotif(vc.getContig(), vc.getStart(), vc.getStart() + size));
    }
}
