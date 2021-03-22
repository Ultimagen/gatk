package org.ultimagenomics.haplotype_calling;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingData;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingLikelihoods;
import org.broadinstitute.hellbender.tools.walkers.genotyper.IndependentSampleGenotypesModel;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.InverseAllele;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.PrintStream;
import java.util.Arrays;

public class AlleleFilteringHC extends AlleleFiltering {
    private HaplotypeCallerGenotypingEngine genotypingEngine;
    public AlleleFilteringHC(HaplotypeCallerArgumentCollection _hcargs, PrintStream assemblyDebugStream, HaplotypeCallerGenotypingEngine _genotypingEngine){
        super(_hcargs, assemblyDebugStream);
        genotypingEngine = _genotypingEngine;
    }


    int getAlleleLikelihood(final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods, Allele allele) {
        final Allele notAllele = InverseAllele.of(allele);

        // iterate over contigs and see what their qual is.

        GenotypingData<Allele> genotypingData = new GenotypingData<>(genotypingEngine.getPloidyModel(), alleleLikelihoods);

        IndependentSampleGenotypesModel genotypesModel = new IndependentSampleGenotypesModel();

        AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(allele, notAllele));

        final GenotypingLikelihoods<Allele> genotypingLikelihoods = genotypesModel.calculateLikelihoods(alleleList, genotypingData);

        final int[] asPL = genotypingLikelihoods.sampleLikelihoods(0).getAsPLs();
        final int retVal;
        retVal = Math.min(asPL[1], asPL[0]) - asPL[2]; // if this is "large", reject the allele.
        logger.debug(() -> String.format("GAL:: %s: %d %d %d", allele.toString(), asPL[0], asPL[1], asPL[2]));
        return retVal;
    }
    
}
