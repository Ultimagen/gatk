package org.broadinstitute.hellbender.tools.walkers.annotator.ultima;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.tools.walkers.annotator.GenotypeAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardMutectAnnotation;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.ultimagen.annotate.AnnotatorBase;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ULTIMA_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ULTIMA_ANNOTATORS_SUMMARY, summary="UltimaClassifyIndel")
public class UltimaClassifyIndel extends UltimaAnnotatorBase {
    private final static Logger logger = LogManager.getLogger(UltimaClassifyIndel.class);

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.ULTIMA_INDEL_CLASSIFY, GATKVCFConstants.ULTIMA_INDEL_LENGTH);
    }


    // "indel_classify" and "indel_length"
    @Override
    public void annotate(final VariantContext vc, final Annotator annotator) {

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
            annotator.annotate(GATKVCFConstants.ULTIMA_INDEL_CLASSIFY, (vc.getReference().length() < maxAlleleLength) ? C_INSERT : C_DELETE);

            /*
            lambda x: max([abs(len(y) - len(x['ref'])) for y in x['alleles']])
             */
            final int refLength = vc.getReference().length();
            annotator.annotate(GATKVCFConstants.ULTIMA_INDEL_LENGTH, vc.getAlleles().stream()
                    .filter(allele -> !allele.isReference())
                    .map(allele -> Math.abs(refLength - allele.length()))
                    .max(Integer::compare).get());
        }
    }

}
