package org.ultimagen.annotate;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public class ClassifyIndelAnnotator extends AnnotatorBase {

    public ClassifyIndelAnnotator(final ReferenceSequenceFile referenceSequenceFile) {
        super(referenceSequenceFile);
    }

    // "indel_classify" and "indel_length"
    @Override
    public void annotate(final VariantContext vc) {

        if ( vc.isIndel() ) {

            /*
            if not x['indel']:
                return None
            elif len(x['ref']) < max([len(y) for y in x['alleles']]):
                return 'ins'
            return 'del'
             */
            final int maxAlleleLength = vc.getAlleles().stream()
                    .filter(allele -> !allele.isReference())
                    .map(allele -> allele.length())
                    .max(Integer::compare).get();
            annotate(vc, A_INDEL_CLASSIFY, (vc.getReference().length() < maxAlleleLength) ? C_INSERT : C_DELETE);

            /*
            lambda x: max([abs(len(y) - len(x['ref'])) for y in x['alleles']])
             */
            final int refLength = vc.getReference().length();
            annotate(vc, A_INDEL_LENGTH, vc.getAlleles().stream()
                    .filter(allele -> !allele.isReference())
                    .map(allele -> Math.abs(refLength - allele.length()))
                    .max(Integer::compare).get());
        }
    }

}
