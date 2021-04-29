package org.ultimagen.haplotypeCalling;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.collections4.ListUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.mutect.*;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.NaturalLogUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
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

        final List<LikelihoodMatrix<GATKRead, Allele>> tumorMatrices = IntStream.range(0, alleleLikelihoods.numberOfSamples())
                .filter(n -> !genotypingEngine.normalSamples.contains(alleleLikelihoods.getSample(n)))
                .mapToObj(alleleLikelihoods::sampleMatrix)
                .collect(Collectors.toList());
        final AlleleList<Allele> alleleList = tumorMatrices.get(0);
        final LikelihoodMatrix<GATKRead, Allele> logTumorMatrix = SomaticGenotypingEngine.combinedLikelihoodMatrix(tumorMatrices, alleleList);
        double normalGermlineLogOdds;
        double normalArtefactLogOdds;
        if (genotypingEngine.normalSamples.size()>0) {
            final List<LikelihoodMatrix<GATKRead, Allele>> normalMatrices = IntStream.range(0, alleleLikelihoods.numberOfSamples())
                    .filter(n -> genotypingEngine.normalSamples.contains(alleleLikelihoods.getSample(n)))
                    .mapToObj(alleleLikelihoods::sampleMatrix)
                    .collect(Collectors.toList());

            final LikelihoodMatrix<GATKRead, Allele> logNormalMatrix = SomaticGenotypingEngine.combinedLikelihoodMatrix(normalMatrices, alleleList);
            normalGermlineLogOdds = diploidAltLogOdds(logNormalMatrix);
            normalArtefactLogOdds = somaticAltLogOdds(logNormalMatrix);
        }else{
            normalGermlineLogOdds = 0;
            normalArtefactLogOdds = 0;
        }

        double tumorLogOdds = somaticAltLogOdds(logTumorMatrix);

        logger.debug(() -> String.format("GAL:: %s: %f %f %f", allele.toString(), tumorLogOdds, normalGermlineLogOdds, normalArtefactLogOdds));
        return (int)(10*tumorLogOdds - 10*Math.max(normalGermlineLogOdds, normalArtefactLogOdds));
    }

    /**
     * Calculate the log likelihoods of the ref/alt het genotype for each alt allele, then subtracts
     * these from the hom ref log likelihood to get the log-odds.
     *
     * @param matrix a matrix of log likelihoods
     */
    private double diploidAltLogOdds(final LikelihoodMatrix<GATKRead, Allele> matrix) {
        final int refIndex = 1;
        final int altIndex = 0;
        final int numReads = matrix.evidenceCount();

        final double homRefLogLikelihood = new IndexRange(0, numReads).sum(r -> matrix.get(refIndex,r));

        final double hetLogLikelihood = new IndexRange(0, numReads)
                .sum(r -> NaturalLogUtils.logSumExp(matrix.get(refIndex, r), matrix.get(altIndex,r)) + NaturalLogUtils.LOG_ONE_HALF);
        final double homAltLogLikelihood = new IndexRange(0, numReads).sum(r -> matrix.get(altIndex,r));
        return homRefLogLikelihood - Math.max(hetLogLikelihood, homAltLogLikelihood);
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
