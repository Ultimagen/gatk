package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.pileup.PileupBasedAlleles;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class HmerIndelBiasAnnotation implements GenotypeAnnotation {

    private final static Logger logger = LogManager.getLogger(HmerIndelBiasAnnotation.class);

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.FLOW_HMER_INDEL_BIAS_KEY);
    }

    @Override
    public void annotate(final ReferenceContext ref, final VariantContext vc, final Genotype g, final GenotypeBuilder gb,
                         final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        // This method should implement the logic to annotate the genotype with Hmer Indel Bias.
        // The actual implementation is not provided in the original code snippet.
        // For now, we will log a message indicating that this method is called.
        List<GATKRead> allReads = likelihoods.sampleEvidence(likelihoods.indexOfSample(g.getSampleName())).stream().collect(Collectors.toList());
        allReads.addAll(likelihoods.filteredSampleEvidence(likelihoods.indexOfSample(g.getSampleName())).stream().collect(Collectors.toList()));
        return;
//        List<GATKRead> leftClippedReads = allReads.stream().filter( rd -> (rd.getStart() <= rd.getEnd()) && rd.overlaps(vc) &&eclowwasReadClipped(rd, false )).collect(Collectors.toList());
//        List<GATKRead> rightClippedReads = allReads.stream().filter( rd -> (rd.getStart() <= rd.getEnd()) && rd.overlaps(vc) && wasReadClipped(rd, true )).collect(Collectors.toList());
//        ReadPileup leftClippedPileup = new ReadPileup(ref.getInterval(),leftClippedReads);
//        Map<Allele, Integer> leftClippedCounts = PileupBasedAlleles.getPileupAlleleCounts(vc, leftClippedPileup);
//        final int[] counts = new int[vc.getNAlleles()];
//        counts[0] = leftClippedCounts.get(vc.getReference()); //first one in AD is always ref
//        for (int i = 0; i < vc.getNAlleles() -1; i++) {
//            counts[i + 1] = leftClippedCounts.get(vc.getAlternateAllele(i));
//        }
//        gb.attribute(getKeyNames().get(0), counts.clone());
//
//        ReadPileup rightClippedPileup = new ReadPileup(ref.getInterval(),rightClippedReads);
//        Map<Allele, Integer> rightClippedCounts = PileupBasedAlleles.getPileupAlleleCounts(vc, rightClippedPileup);
//        for (int i = 0 ; i < counts.length; i++){
//            counts[i] = 0;
//        }
//        counts[0] = rightClippedCounts.get(vc.getReference()); //first one in AD is always ref
//        for (int i = 0; i < vc.getNAlleles() -1; i++) {
//            counts[i + 1] = rightClippedCounts.get(vc.getAlternateAllele(i));
//        }
//
    }
}
