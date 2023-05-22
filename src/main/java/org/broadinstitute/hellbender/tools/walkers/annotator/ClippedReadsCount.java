package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AlleleSpecificAnnotationData;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

/**
 * Counts number of soft-clipped reads that support each allele. The format is
 * #left clips, #right clips | #left clips, #right clips | ... for each allele
 *
 * <h3> Caveats </h3>
 * This annotation can only be calculated by Mutect2 and HaplotypeCaller. In case
 * the soft clipping is reverted during the flow
 */

@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Number of s")

public class ClippedReadsCount implements GenotypeAnnotation {
    private final static Logger logger = LogManager.getLogger(ClippedReadsCount.class);

    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.SCL_COUNT_KEY, GATKVCFConstants.SCR_COUNT_KEY);
    }

    @Override
    public VCFCompoundHeaderLine.SupportedHeaderLineType annotationType() {
        return GenotypeAnnotation.super.annotationType();
    }

    @Override
    public List<VCFCompoundHeaderLine> getDescriptions() {
        return GenotypeAnnotation.super.getDescriptions();
    }

    @Override
    public void annotate(ReferenceContext ref, VariantContext vc, Genotype g, GenotypeBuilder gb, AlleleLikelihoods<GATKRead, Allele> likelihoods) {

    }
}
