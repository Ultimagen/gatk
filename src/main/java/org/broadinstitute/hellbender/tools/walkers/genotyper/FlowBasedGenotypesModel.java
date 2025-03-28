package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.utils.dragstr.DragstrReferenceAnalyzer;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.AlleleListPermutation;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.List;

/**
 * This is the genotyper model that genotypes homopolymer length.
 *
 *
 */

public class FlowBasedGenotypesModel extends DRAGENGenotypesModel {

    final double refBias;
    public FlowBasedGenotypesModel(final boolean useBQDModel, final boolean useFRDModel, final int allelePadding,
                                   final int maxEffectiveDepthAdjustment,
                                   final double maxForeignReadFraction,
                                   final DragstrParams dragstrParams,
                                   final double refBias) {
        super(useBQDModel, useFRDModel, allelePadding,
        maxEffectiveDepthAdjustment,
        maxForeignReadFraction,
        dragstrParams);
        this.refBias = refBias;
    }

    @Override
    public <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(AlleleList<A> genotypingAlleles, GenotypingData<A> data,
                                                                            byte[] paddedReference, int offsetForRefIntoEvent,
                                                                            DragstrReferenceAnalyzer dragstrs) {
        GenotypingLikelihoods<A> dragenLikelihoods = super.calculateLikelihoods(genotypingAlleles, data, paddedReference, offsetForRefIntoEvent, dragstrs);

        return dragenLikelihoods;
    }

    public <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(final AlleleList<A> genotypingAlleles,
                                                                            final GenotypingData<A> data,
                                                                            final byte[] paddedReference,
                                                                            final int offsetForRefIntoEvent,
                                                                            final DragstrReferenceAnalyzer dragstrs,
                                                                            final boolean isRefBiasExpected) {
        Utils.nonNull(genotypingAlleles, "the allele cannot be null");
        Utils.nonNull(data, "the genotyping data cannot be null");

        final AlleleListPermutation<A> permutation = data.permutation(genotypingAlleles);
        final AlleleLikelihoodMatrixMapper<A> alleleLikelihoodMatrixMapper = new AlleleLikelihoodMatrixMapper<>(permutation);

        final int sampleCount = data.numberOfSamples();
        final PloidyModel ploidyModel = data.ploidyModel();
        GenotypingLikelihoods<A> genotypeLikelihoods = null;
        final List<GenotypeLikelihoods> likelihoods = new ArrayList<>(sampleCount);

        if (data.readLikelihoods().getVariantCallingSubsetApplied()!=null) { // in the context of AlleleFiltering
            genotypeLikelihoods = super.calculateLikelihoods(genotypingAlleles, data, paddedReference, offsetForRefIntoEvent, dragstrs);
        } else {
            for (int i = 0; i < sampleCount; i++) {
                final int samplePloidy = ploidyModel.samplePloidy(i);
                final LikelihoodMatrix<GATKRead, A> sampleLikelihoods = alleleLikelihoodMatrixMapper.mapAlleles(data.readLikelihoods().sampleMatrix(i));
                likelihoods.add(GenotypeLikelihoodCalculator.log10GenotypeLikelihoods(samplePloidy, sampleLikelihoods, refBias));
            }
            genotypeLikelihoods = new GenotypingLikelihoods<>(genotypingAlleles, ploidyModel, likelihoods);
        }
        for (int i = 0; i < sampleCount; i++) {
            final int samplePloidy = ploidyModel.samplePloidy(i);

            final LikelihoodMatrix<GATKRead, A> sampleLikelihoods = alleleLikelihoodMatrixMapper.mapAlleles(data.readLikelihoods().sampleMatrix(i));
            if (isRefBiasExpected) {
                applyLikelihoodsAdjusmentToBaseline(genotypeLikelihoods.sampleLikelihoods(i).getAsVector(), "BIAS",
                        GenotypeLikelihoodCalculator.log10GenotypeLikelihoods(samplePloidy, sampleLikelihoods, refBias).getAsVector());
            }
        }

        return genotypeLikelihoods;
    }

    static public boolean isEligibleHomopolymerIndel(final VariantContext vc, final int loc,
                                                     final DragstrReferenceAnalyzer dragstrs, final int hpolIndelThreshold) {
        if (vc==null) {
            return false;
        }
        final int period = dragstrs.period(loc);
        final int repeats = dragstrs.repeatLength(loc);
        final byte ru = dragstrs.repeatUnit(loc)[0];
        if ((period == 1) && (repeats >= hpolIndelThreshold)){
            if (!vc.isIndel() || !vc.getAlleles().stream().allMatch(a -> isHmerIndel(a,ru))) {
                return false;
            }
        } else {
            return false;
        }
        return true;
    }

    static public boolean isEligibleHomopolymerIndel(final Event allele, final int loc,
                                                     final DragstrReferenceAnalyzer dragstrs, final int hpolIndelThreshold) {
        if (allele==null) {
            return false;
        }
        final int period = dragstrs.period(loc);
        final int repeats = dragstrs.repeatLength(loc);
        final byte ru = dragstrs.repeatUnit(loc)[0];
        if ((period == 1) && (repeats >= hpolIndelThreshold)){
            if (!allele.isIndel() || !isHmerIndel(allele.refAllele(), ru) || !isHmerIndel(allele.altAllele(), ru)) {
                return false;
            }
        } else {
            return false;
        }
        return true;
    }


    static protected boolean isHmerIndel(final Allele al, final byte hmer_base){
        for (int i = 1; i< al.length(); i++){
            if (al.getBases()[i] != hmer_base){
                return false;
            }
        }
        return true;
    }
    public static boolean isRefBiasExpected(final Event allele,
                                            final DragstrReferenceAnalyzer dragstrReferenceAnalyzer,
                                            final int alleleOffset,
                                            final int biasedIndelThreshold) {
        if (biasedIndelThreshold == 0) {
            return false;
        }
        if (!allele.isSimpleInsertion()) {
            return false;
        }
        return isEligibleHomopolymerIndel(allele, alleleOffset, dragstrReferenceAnalyzer, biasedIndelThreshold);
    }
}
