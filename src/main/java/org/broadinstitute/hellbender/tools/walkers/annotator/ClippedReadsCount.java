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
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceModel;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Counts number of soft-clipped reads that support each allele.  *
 * <h3> Caveats </h3>
 * This annotation can only be calculated by Mutect2 and HaplotypeCaller. In case
 * the soft clipping is reverted during the flow
 */

@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Number of soft clipped reads ")

public class ClippedReadsCount implements GenotypeAnnotation {
    private final static Logger logger = LogManager.getLogger(ClippedReadsCount.class);

    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.SOFT_CLIP_LEFT_COUNT_KEY, GATKVCFConstants.SOFT_CLIP_RIGHT_COUNT_KEY);
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

        List<GATKRead> allReads = likelihoods.sampleEvidence(likelihoods.indexOfSample(g.getSampleName())).stream().collect(Collectors.toList());

        List<GATKRead> leftClipped = allReads.stream().filter( rd -> rd.hasAttribute(ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_START_TAG)).collect(Collectors.toList());
        List<GATKRead> rightClipped = allReads.stream().filter( rd -> rd.hasAttribute(ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_END_TAG)).collect(Collectors.toList());
    }
}
