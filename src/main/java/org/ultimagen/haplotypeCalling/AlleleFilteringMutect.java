package org.ultimagen.haplotypeCalling;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.walkers.mutect.*;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.OutputStreamWriter;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class AlleleFilteringMutect extends AlleleFiltering {
    private SomaticGenotypingEngine genotypingEngine;
    public AlleleFilteringMutect(M2ArgumentCollection _m2args,
                                 OutputStreamWriter assemblyDebugStream,
                                 SomaticGenotypingEngine _genotypingEngine){
        super(_m2args, assemblyDebugStream);
        genotypingEngine = _genotypingEngine;
    }


    int getAlleleLikelihood(final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods, Allele allele) {

        final List<LikelihoodMatrix<GATKRead, Allele>> allMatrices = IntStream.range(0, alleleLikelihoods.numberOfSamples())
                .mapToObj(alleleLikelihoods::sampleMatrix)
                .collect(Collectors.toList());
        final AlleleList<Allele> alleleList = allMatrices.get(0);
        final LikelihoodMatrix<GATKRead, Allele> logAllMatrix = SomaticGenotypingEngine.combinedLikelihoodMatrix(allMatrices, alleleList);
        double alleleLogOdds = somaticAltLogOdds(logAllMatrix);
        logger.debug(() -> String.format("GAL:: %s: %f", allele.toString(), alleleLogOdds));
        return (int)(10*alleleLogOdds);
    }


    /**
     * Calculate the log likelihoods of the ref/alt het genotype for each alt allele, then subtracts
     * these from the hom ref log likelihood to get the log-odds.
     *
     * @param matrix a matrix of log likelihoods
     */
    private double somaticAltLogOdds(final LikelihoodMatrix<GATKRead, Allele> matrix) {

        final LikelihoodMatrix<GATKRead, Allele> initialMatrix = matrix;
        final LikelihoodMatrix<GATKRead, Allele> logMatrixWithoutThisAllele = SubsettedLikelihoodMatrix.excludingAllele(matrix, matrix.getAllele(0));

        final double logEvidenceWithoutThisAllele = logMatrixWithoutThisAllele.evidenceCount() == 0 ? 0 :
                SomaticLikelihoodsEngine.logEvidence(SomaticGenotypingEngine.getAsRealMatrix(logMatrixWithoutThisAllele),
                        genotypingEngine.makePriorPseudocounts(logMatrixWithoutThisAllele.numberOfAlleles()));
        final double logEvidenceWithAllAlleles= initialMatrix.evidenceCount() == 0 ? 0 :
                SomaticLikelihoodsEngine.logEvidence(SomaticGenotypingEngine.getAsRealMatrix(initialMatrix),
                        genotypingEngine.makePriorPseudocounts(initialMatrix.numberOfAlleles()));
        double tumorLogOdds = (-logEvidenceWithAllAlleles + logEvidenceWithoutThisAllele);
        return tumorLogOdds;
    }
}
