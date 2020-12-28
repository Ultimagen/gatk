package org.ultimagenomics.haplotype_calling;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandOddsRatio;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingData;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingLikelihoods;
import org.broadinstitute.hellbender.tools.walkers.genotyper.IndependentSampleGenotypesModel;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.InverseAllele;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.JoinedContigs;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.PrintStream;
import java.util.Arrays;

public class ContigFilteringHC extends ContigFiltering {
    private HaplotypeCallerGenotypingEngine genotypingEngine;
    public ContigFilteringHC(HaplotypeCallerArgumentCollection _hcargs, PrintStream assemblyDebugStream, HaplotypeCallerGenotypingEngine _genotypingEngine){
        super(_hcargs, assemblyDebugStream);
        genotypingEngine = _genotypingEngine;
    }


    int getContigLikelihood(final AlleleLikelihoods<GATKRead, Allele> contigLikelihoods, Allele contig) {
        final Allele notContig = InverseAllele.of(contig);

        // iterate over contigs and see what their qual is.

        GenotypingData<Allele> genotypingData = new GenotypingData<>(genotypingEngine.getPloidyModel(), contigLikelihoods);

        IndependentSampleGenotypesModel genotypesModel = new IndependentSampleGenotypesModel();

        AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(contig, notContig));

        final GenotypingLikelihoods<Allele> genotypingLikelihoods = genotypesModel.calculateLikelihoods(alleleList, genotypingData);

        final int[] asPL = genotypingLikelihoods.sampleLikelihoods(0).getAsPLs();
        final int retVal;
        retVal = Math.min(asPL[1], asPL[0]) - asPL[2]; // if this is "large", reject the contig.
        logger.debug(String.format("GCL:: %s->%s: %d %d %d", ((JoinedContigs)contig).getAllele1().getBaseString(),
                ((JoinedContigs)contig).getAllele2().getBaseString(), asPL[0], asPL[1], asPL[2]));
        return retVal;
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
