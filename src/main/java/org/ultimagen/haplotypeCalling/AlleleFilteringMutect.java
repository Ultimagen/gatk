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
        final LikelihoodMatrix<GATKRead, Allele> logTumorMatrix = genotypingEngine.combinedLikelihoodMatrix(tumorMatrices, alleleList);

        final List<LikelihoodMatrix<GATKRead, Allele>> normalMatrices = IntStream.range(0, alleleLikelihoods.numberOfSamples())
                .filter(n -> genotypingEngine.normalSamples.contains(alleleLikelihoods.getSample(n)))
                .mapToObj(alleleLikelihoods::sampleMatrix)
                .collect(Collectors.toList());
        final LikelihoodMatrix<GATKRead, Allele> logNormalMatrix = genotypingEngine.combinedLikelihoodMatrix(normalMatrices, alleleList);


        final LikelihoodMatrix<GATKRead, Allele> initialMatrix = logTumorMatrix;
        final LikelihoodMatrix<GATKRead, Allele> logMatrixWithoutThisAllele = SubsettedLikelihoodMatrix.excludingAllele(logTumorMatrix, allele);
        final double normalLogOdds = diploidAltLogOdds(logNormalMatrix);


        final double logEvidenceWithoutThisAllele = logMatrixWithoutThisAllele.evidenceCount() == 0 ? 0 :
                SomaticLikelihoodsEngine.logEvidence(SomaticGenotypingEngine.getAsRealMatrix(logMatrixWithoutThisAllele),
                        genotypingEngine.makePriorPseudocounts(logMatrixWithoutThisAllele.numberOfAlleles()));
        final double logEvidenceWithAllAlleles= initialMatrix.evidenceCount() == 0 ? 0 :
                SomaticLikelihoodsEngine.logEvidence(SomaticGenotypingEngine.getAsRealMatrix(initialMatrix),
                        genotypingEngine.makePriorPseudocounts(initialMatrix.numberOfAlleles()));

        logger.debug(() -> String.format("GAL:: %s: %f %f %d", allele.toString(), logEvidenceWithAllAlleles,
                logEvidenceWithoutThisAllele,
                (int)(10*(-logEvidenceWithAllAlleles + logEvidenceWithoutThisAllele))));

        return (int)(10*(-logEvidenceWithAllAlleles + logEvidenceWithoutThisAllele));



        final PerAlleleCollection<Double> tumorLogOdds = somaticLogOdds(logTumorMatrix);

        final PerAlleleCollection<Double> normalArtifactLogOdds = somaticLogOdds(logNormalMatrix);


        final Set<Allele> forcedAlleles = AssemblyBasedCallerUtils.getAllelesConsistentWithGivenAlleles(givenAlleles, mergedVC);
        final List<Allele> tumorAltAlleles = mergedVC.getAlternateAlleles().stream()
                .filter(allele -> forcedAlleles.contains(allele) || tumorLogOdds.getAlt(allele) > MTAC.getEmissionLogOdds())
                .collect(Collectors.toList());

        final long somaticAltCount = tumorAltAlleles.stream()
                .filter(allele -> forcedAlleles.contains(allele) || !hasNormal || MTAC.genotypeGermlineSites || normalLogOdds.getAlt(allele) > MathUtils.log10ToLog(MTAC.normalLog10Odds))
                .count();

        // if every alt allele is germline, skip this variant.  However, if some alt alleles are germline and others
        // are not we emit them all so that the filtering engine can see them
        if (somaticAltCount == 0) {
            continue;
        }

        final List<Allele> allAllelesToEmit = ListUtils.union(Arrays.asList(mergedVC.getReference()), tumorAltAlleles);




    }

    /**
     * Calculate the log likelihoods of the ref/alt het genotype for each alt allele, then subtracts
     * these from the hom ref log likelihood to get the log-odds.
     *
     * @param matrix a matrix of log likelihoods
     */
    private double diploidAltLogOdds(final LikelihoodMatrix<GATKRead, Allele> matrix) {
        final int refIndex = 0;
        final int altIndex = 1;
        final int numReads = matrix.evidenceCount();

        final double homRefLogLikelihood = new IndexRange(0, numReads).sum(r -> matrix.get(refIndex,r));

        final double hetLogLikelihood = new IndexRange(0, numReads)
                .sum(r -> NaturalLogUtils.logSumExp(matrix.get(refIndex, r), matrix.get(altIndex,r)) + NaturalLogUtils.LOG_ONE_HALF);
        final double homAltLogLikelihood = new IndexRange(0, numReads).sum(r -> matrix.get(altIndex,r));
        return homRefLogLikelihood - Math.max(hetLogLikelihood, homAltLogLikelihood);
    }

}
