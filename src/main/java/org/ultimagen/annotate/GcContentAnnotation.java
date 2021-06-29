package org.ultimagen.annotate;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public class GcContentAnnotation extends AnnotatorBase {

    private int windowSize = 10;

    public GcContentAnnotation(final ReferenceSequenceFile referenceSequenceFile) {
        super(referenceSequenceFile);
    }

    @Override
    public void annotate(final VariantContext vc) {

        /*
        chrom = faidx[rec['chrom']]
        beg = rec['pos'] - int(size / 2)
        end = beg + size
        seq = chrom[beg:end].seq.upper()
        seqGC = seq.replace('A', '').replace('T', '')
        return float(len(seqGC)) / len(seq)
         */
        int         begin = vc.getStart() - (int)(windowSize / 2);
        int         end = begin + windowSize;
        String      seq = getReferenceMotif(vc.getContig(), begin, end);
        int         gcCount = 0;
        for ( byte b : seq.getBytes() ) {
            if ( b == (char)'G' || b == (char)'C' ) {
                gcCount++;
            }
        }
        annotate(vc, A_GC_CONTENT, (float)gcCount / seq.length());
    }
}
