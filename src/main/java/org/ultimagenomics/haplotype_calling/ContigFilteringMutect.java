package org.ultimagenomics.haplotype_calling;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandOddsRatio;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.InverseAllele;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.JoinedContigs;
import org.broadinstitute.hellbender.tools.walkers.mutect.M2ArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.SomaticGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.SomaticLikelihoodsEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.SubsettedLikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.PrintStream;
import java.util.Arrays;


public class ContigFilteringMutect extends ContigFiltering {
    private SomaticGenotypingEngine genotypingEngine;
    public ContigFilteringMutect(M2ArgumentCollection _m2args, PrintStream assemblyDebugStream, SomaticGenotypingEngine _genotypingEngine){
        super(_m2args, assemblyDebugStream);
        genotypingEngine = _genotypingEngine;
    }

    int getContigLikelihood(final AlleleLikelihoods<GATKRead, Allele> contigLikelihoods, Allele contig) {

        final LikelihoodMatrix<GATKRead, Allele> logMatrixWithoutThisAllele = SubsettedLikelihoodMatrix.excludingAllele(contigLikelihoods.sampleMatrix(0), contig);
        final double logEvidenceWithoutThisAllele = logMatrixWithoutThisAllele.evidenceCount() == 0 ? 0 :
                SomaticLikelihoodsEngine.logEvidence(SomaticGenotypingEngine.getAsRealMatrix(logMatrixWithoutThisAllele),
                        genotypingEngine.makePriorPseudocounts(logMatrixWithoutThisAllele.numberOfAlleles()));

        return (int)(10*logEvidenceWithoutThisAllele);

    }

    double getContigSOR(final AlleleLikelihoods<GATKRead, Allele> contigLikelihoods, Allele contig) {
        final Allele notContig = InverseAllele.of(contig);
        int [][] contingency_table = StrandOddsRatio.getContingencyTable(contigLikelihoods, notContig, Arrays.asList(contig), 1);
        double sor = StrandOddsRatio.calculateSOR(contingency_table);
        logger.debug(String.format("GCS:: %s->%s: %f (%d %d %d %d)", ((JoinedContigs)contig).getAllele1().getBaseString(),
                ((JoinedContigs)contig).getAllele2().getBaseString(), sor, contingency_table[0][0], contingency_table[0][1], contingency_table[1][0], contingency_table[1][1]));
        return sor;

    }
}
