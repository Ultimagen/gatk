package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.apache.commons.lang3.NotImplementedException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrReferenceAnalyzer;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.AlleleListPermutation;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.List;

public class IndependentSampleGenotypesModelWithRefBias implements GenotypingModel{
    public IndependentSampleGenotypesModelWithRefBias() { }

    public <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(final AlleleList<A> genotypingAlleles,
                                                                            final GenotypingData<A> data,
                                                                            final byte[] paddedReference,
                                                                            final int offsetForRefIntoEvent,
                                                                            final DragstrReferenceAnalyzer dragstrs) throws NotImplementedException {
        throw new NotImplementedException("calculateLikelihoods without reference bias not implemented for this model");
    }
    public <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(final AlleleList<A> genotypingAlleles,
                                                                            final GenotypingData<A> data,
                                                                            final byte[] paddedReference,
                                                                            final int offsetForRefIntoEvent,
                                                                            final double refBias) {
        Utils.nonNull(genotypingAlleles, "the allele cannot be null");
        Utils.nonNull(data, "the genotyping data cannot be null");

        final AlleleListPermutation<A> permutation = data.permutation(genotypingAlleles);
        final AlleleLikelihoodMatrixMapper<A> alleleLikelihoodMatrixMapper = new AlleleLikelihoodMatrixMapper<>(permutation);

        final int sampleCount = data.numberOfSamples();
        final PloidyModel ploidyModel = data.ploidyModel();
        final List<GenotypeLikelihoods> genotypeLikelihoods = new ArrayList<>(sampleCount);

        for (int i = 0; i < sampleCount; i++) {
            final int samplePloidy = ploidyModel.samplePloidy(i);

            final LikelihoodMatrix<GATKRead, A> sampleLikelihoods = alleleLikelihoodMatrixMapper.mapAlleles(data.readLikelihoods().sampleMatrix(i));
            genotypeLikelihoods.add(GenotypeLikelihoodCalculator.log10GenotypeLikelihoods(samplePloidy, sampleLikelihoods, refBias));
        }
        return new GenotypingLikelihoods<>(genotypingAlleles, ploidyModel, genotypeLikelihoods);
    }
}

