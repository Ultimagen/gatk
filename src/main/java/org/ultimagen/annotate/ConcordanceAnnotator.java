package org.ultimagen.annotate;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;

public class ConcordanceAnnotator implements VariantContextAnnotator {

    private List<VariantContextAnnotator>   annotators = new LinkedList<>();

    public ConcordanceAnnotator(final ReferenceSequenceFile referenceSequenceFile) {
        annotators.add(new ClassifyIndelAnnotator(referenceSequenceFile));
        annotators.add(new IsHmerIndelAnnotator(referenceSequenceFile));
        annotators.add(new GetMotifAroundAnnotation(referenceSequenceFile));
    }

    @Override
    public void annotate(final VariantContext vc) {
        annotators.forEach(a -> a.annotate(vc));
    }
}
