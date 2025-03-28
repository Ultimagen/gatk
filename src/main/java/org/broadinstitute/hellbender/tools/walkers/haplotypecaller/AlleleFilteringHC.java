package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculationResult;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AlleleFrequencyCalculator;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.InverseAllele;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Filtering haplotypes that contribute weak alleles to the genotyping. This is a version that determines if allele is weak using
 * HaplotypeCaller model.
 *
 * @author Ilya Soifer &lt;ilya.soifer@ultimagen.com&gt;
 * @author Yossi Farjoun &lt;farjoun@broadinstitute.org&gt;
 *
 */

public class AlleleFilteringHC extends AlleleFiltering {
    private HaplotypeCallerGenotypingEngine genotypingEngine;
    private AlleleFrequencyCalculator afCalc;
    final double insertionRefBias;
    public AlleleFilteringHC(HaplotypeCallerArgumentCollection _hcargs, OutputStreamWriter assemblyDebugStream,
                             HaplotypeCallerGenotypingEngine _genotypingEngine, final SAMFileHeader header) {
        super(_hcargs, assemblyDebugStream, header, _hcargs == null ? 0:_hcargs.homopolymerGenotypingThreshold);
        genotypingEngine = _genotypingEngine;
        GenotypeCalculationArgumentCollection config = genotypingEngine.getConfiguration().genotypeArgs;
         afCalc = AlleleFrequencyCalculator.makeCalculator(config);
         if (_hcargs==null){
             this.insertionRefBias = 0.5;
         } else {
             this.insertionRefBias = _hcargs.insertionRefBias;
         }
    }

    protected double getStringentQuality() { return 1; }
    /**
     * Calculate genotype likelihood of requirement of an allele. Specifically, calculates the likelihood
     * of the data given that allele versus the likelihood of the data when all haplotypes containing the allele are removed
     * This is very similar to what is done in the genotyping engine, but here the haplotypes that do not support the allele
     * are all haplotypes that do not contain the allele.
     *
     * @param alleleLikelihoods
     * @param allele
     * @return likelihood, expressed as phred-scaled confidence
     */
    @Override
    int getAlleleLikelihoodVsInverse(final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods, final Allele allele, final boolean isRefBiasExpected) {
        final Allele notAllele = InverseAllele.of(allele, true);

        // iterate over contigs and see what their qual is.

        GenotypingData<Allele> genotypingData = new GenotypingData<>(genotypingEngine.getPloidyModel(), alleleLikelihoods);

        final FlowBasedGenotypesModel genotypesModel = new FlowBasedGenotypesModel(false, false,
                assemblyArgs.informativeReadOverlapMargin, 0,
                1, null, insertionRefBias );

        AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(notAllele, allele));

        final GenotypingLikelihoods<Allele> genotypingLikelihoods = genotypesModel.calculateLikelihoods(alleleList,
                genotypingData, null, 0, null, isRefBiasExpected);

        List<Integer> perSamplePLs = new ArrayList<>();
        for (int i = 0; i < genotypingLikelihoods.numberOfSamples(); i++) {
            final int[] pls = genotypingLikelihoods.sampleLikelihoods(i).getAsPLs();
            perSamplePLs.add(Math.min(pls[1] - pls[0], pls[2] - pls[0]));

            final int finalI = i;
            logger.debug(() -> String.format("GAL (%s):: %s: %d %d %d",
                    genotypingLikelihoods.getSample(finalI), allele.toString(), pls[0], pls[1], pls[2]));
        }
        return Collections.min(perSamplePLs);
    }
}
