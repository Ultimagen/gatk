package org.ultimagen.haplotypeCalling;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.SomaticGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.SomaticLikelihoodsEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.SubsettedLikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.PrintStream;


public class AlleleFilteringMutect extends AlleleFiltering {
    private SomaticGenotypingEngine genotypingEngine;
    public AlleleFilteringMutect(M2ArgumentCollection _m2args, PrintStream assemblyDebugStream, SomaticGenotypingEngine _genotypingEngine){
        super(_m2args, assemblyDebugStream);
        genotypingEngine = _genotypingEngine;
    }

    int getAlleleLikelihood(final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods, Allele allele) {
        final LikelihoodMatrix<GATKRead, Allele> initialMatrix = alleleLikelihoods.sampleMatrix(0);
        final LikelihoodMatrix<GATKRead, Allele> logMatrixWithoutThisAllele = SubsettedLikelihoodMatrix.excludingAllele(alleleLikelihoods.sampleMatrix(0), allele);


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

    }

}
