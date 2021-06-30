package org.broadinstitute.hellbender.tools.walkers.annotator.ultima;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.GenotypeAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardMutectAnnotation;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.Collections;
import java.util.List;

@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ULTIMA_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ULTIMA_ANNOTATORS_SUMMARY, summary="Development annotation 1")
public final class DevelAnnotation1 extends GenotypeAnnotation implements StandardMutectAnnotation {
    private final static Logger logger = LogManager.getLogger(DevelAnnotation1.class);

    @Override
    public void annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        Utils.nonNull(vc);
        Utils.nonNull(g);
        Utils.nonNull(gb);

        if ( likelihoods == null || !g.isCalled() ) {
            logger.warn(AnnotationUtils.generateMissingDataWarning(vc, g, likelihoods));
            return;
        }

        gb.attribute(GATKVCFConstants.ULTIMA_DEVEL_ANNOTATION_1, 1111);
    }

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GATKVCFConstants.ULTIMA_DEVEL_ANNOTATION_1);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Collections.singletonList(GATKVCFHeaderLines.getFormatLine(getKeyNames().get(0)));
    }

}
