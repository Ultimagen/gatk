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
import org.broadinstitute.hellbender.tools.funcotator.FilterFuncotations;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.pileup.PileupBasedAlleles;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;

@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Mean counts of mismatches per allele")
public class MismatchCount implements GenotypeAnnotation {
    private final static Logger logger = LogManager.getLogger(MismatchCount.class);

    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.NM_COUNT_KEY);
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

        Map<Allele, List<Integer>> mismatchCounts = new HashMap<>();
        for (Allele al: vc.getAlleles()){
            mismatchCounts.put(al, new ArrayList<>());
        }

        fillMDFromLikelihoods(vc, ref, likelihoods, mismatchCounts);
        final int[] counts = new int[vc.getNAlleles()];

        if (mismatchCounts.get(vc.getReference()).size()==0){
            counts[0] = 0;
        } else {
            counts[0] = (int) MathUtils.median(mismatchCounts.get(vc.getReference()));
        }
        for (int i = 0; i < vc.getNAlleles() -1; i++) {
            if (mismatchCounts.get(vc.getAlternateAllele(i)).size()==0){
                counts[i+1] = 0;
            } else {
                counts[i + 1] = (int) MathUtils.median(mismatchCounts.get(vc.getAlternateAllele(i)));
            }
        }
        gb.attribute(getKeyNames().get(0), counts);
    }


    protected void fillMDFromLikelihoods(VariantContext vc, ReferenceContext ref, AlleleLikelihoods<GATKRead, Allele> likelihoods, Map<Allele,List<Integer>> mismatchCounts ) {
        for (final AlleleLikelihoods<GATKRead, Allele>.BestAllele bestAllele : likelihoods.bestAllelesBreakingTies()) {
            final GATKRead read = bestAllele.evidence;
            final Allele allele = bestAllele.allele;
            if (bestAllele.isInformative() && isUsableRead(read, vc)) {
                final Integer value = getElementForRead(read, vc, ref);
                mismatchCounts.get(allele).add(value);
            }
        }
    }
    protected Integer getElementForRead(final GATKRead read, final VariantContext vc, final ReferenceContext ref){
        return AlignmentUtils.getMismatchCount(read, ref.getBases(), read.getStart()).numMismatches;
    }
    /**
     * Can the read be used in comparative tests between ref / alt bases?
     *
     * @param read   the read to consider
     * @param vc    the variant to be annotated
     * @return true if this read is meaningful for comparison, false otherwise
     */
    protected boolean isUsableRead(final GATKRead read, final VariantContext vc) {
        Utils.nonNull(read);
        return read.getMappingQuality() != 0 && read.getMappingQuality() != QualityUtils.MAPPING_QUALITY_UNAVAILABLE;
    }

}