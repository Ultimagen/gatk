package org.ultimagen.haplotypeCalling;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import org.apache.commons.math3.util.MathArrays;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculationResult;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AlleleFrequencyCalculator;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.InverseAllele;
import org.broadinstitute.hellbender.utils.Dirichlet;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class AlleleFilteringHC extends AlleleFiltering {
    private HaplotypeCallerGenotypingEngine genotypingEngine;
    private AlleleFrequencyCalculator afCalc;

    public AlleleFilteringHC(HaplotypeCallerArgumentCollection _hcargs, OutputStreamWriter assemblyDebugStream, HaplotypeCallerGenotypingEngine _genotypingEngine){
        super(_hcargs, assemblyDebugStream);
        genotypingEngine = _genotypingEngine;
        GenotypeCalculationArgumentCollection config = genotypingEngine.getConfiguration().genotypeArgs;
         afCalc = AlleleFrequencyCalculator.makeCalculator(config);
    }


    int getAlleleLikelihood(final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods, Allele allele) {
        final Allele notAllele = InverseAllele.of(allele, true);

        // iterate over contigs and see what their qual is.

        GenotypingData<Allele> genotypingData = new GenotypingData<>(genotypingEngine.getPloidyModel(), alleleLikelihoods);

        IndependentSampleGenotypesModel genotypesModel = new IndependentSampleGenotypesModel();

        AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(notAllele, allele));

        final GenotypingLikelihoods<Allele> genotypingLikelihoods = genotypesModel.calculateLikelihoods(alleleList,
                genotypingData, null, 0, null);
        AFCalculationResult af = afCalc.calculate(genotypingLikelihoods, genotypingEngine.getPloidyModel().totalPloidy());
        final double log10Confidence = af.log10ProbOnlyRefAlleleExists();
        final double phredScaledConfidence = (10.0 * log10Confidence) + 0.0;

        final int[] asPL = genotypingLikelihoods.sampleLikelihoods(0).getAsPLs();

        logger.debug(() -> String.format("GAL:: %s: %d %d %d", allele.toString(), asPL[0], asPL[1], asPL[2]));
        return (int)phredScaledConfidence;
    }

}
