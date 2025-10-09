package org.broadinstitute.hellbender.tools.walkers.annotator.flow;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.GenotypeAnnotation;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.pileup.FlowBasedPileupElement;
import org.broadinstitute.hellbender.utils.pileup.FlowBasedReadPileup;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;

public class HmerIndelBias extends FlowAnnotatorBase implements GenotypeAnnotation {

    private static final OneShotLogger oneShotLogger = new OneShotLogger(HmerIndelBias.class);

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.FLOW_HMER_INDEL_BIAS_KEY);
    }

    @Override
    public htsjdk.variant.vcf.VCFCompoundHeaderLine.SupportedHeaderLineType annotationType() {
        return htsjdk.variant.vcf.VCFCompoundHeaderLine.SupportedHeaderLineType.FORMAT;
    }

    /**
     * Determines if the input VariantContext is a homopolymer indel
     * @param vc the VariantContext to check
     * @param ref the reference context
     * @param likelihoods the allele likelihoods
     * @return true if the variant is a homopolymer indel, false otherwise
     */
    private boolean isHomopolymerIndel(VariantContext vc, ReferenceContext ref, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        // Create LocalContext using FlowAnnotatorBase's infrastructure
        final LocalContext localContext = new LocalContext(ref, vc, likelihoods, true);
        
        if (!localContext.generateAnnotation) {
            return false;
        }
        
        // Use existing FlowAnnotatorBase methods to classify the variant
        indelClassify(vc, localContext);
        isHmerIndel(vc, localContext);
        variantType(vc, localContext);
        
        // Check if variant type is homopolymer indel
        Object variantType = localContext.attributes.get(GATKVCFConstants.FLOW_VARIANT_TYPE);
        return C_H_MER.equals(variantType);
    }

    // Required by InfoFieldAnnotation interface (inherited through FlowAnnotatorBase)
    @Override
    public Map<String, Object> annotate(final ReferenceContext ref,
                                        final VariantContext vc,
                                        final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        // This method is required by InfoFieldAnnotation but not used for GenotypeAnnotation
        return java.util.Collections.emptyMap();
    }

    @Override
    public void annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        // If the input VariantContext is not homopolymer indel - do nothing
        if (!isHomopolymerIndel(vc, ref, likelihoods)) {
            return;
        }
        
        // Get the sample name from the genotype
        final String sampleName = g.getSampleName();
        
        // Find the sample index in the likelihoods
        int sampleIndex = -1;
        for (int i = 0; i < likelihoods.numberOfSamples(); i++) {
            if (likelihoods.getSample(i).equals(sampleName)) {
                sampleIndex = i;
                break;
            }
        }
        
        if (sampleIndex == -1) {
            oneShotLogger.warn("Sample " + sampleName + " not found in likelihoods for variant at " + vc.getContig() + ":" + vc.getStart());
            return;
        }
        
        // Collect all reads (both evidence and filtered evidence) for this sample
        List<GATKRead> allReads = likelihoods.sampleEvidence(sampleIndex).stream().collect(Collectors.toList());
        allReads.addAll(likelihoods.filteredSampleEvidence(sampleIndex).stream().collect(Collectors.toList()));

        // Filter for FlowBasedReads only
        List<FlowBasedRead> flowBasedReads = allReads.stream()
            .filter(read -> read instanceof FlowBasedRead)
            .map(read -> (FlowBasedRead) read)
            .collect(Collectors.toList());

        // Check if we have any flow-based reads
        if (flowBasedReads.isEmpty()) {
            oneShotLogger.warn("No flow-based reads found for sample " + sampleName + " at variant " + vc.getContig() + ":" + vc.getStart());
            return;
        }

        // Create FlowBasedReadPileup
        FlowBasedReadPileup flowPileup = new FlowBasedReadPileup(ref.getInterval(), flowBasedReads);
        
        // Get homopolymer information using FlowAnnotatorBase
        final LocalContext localContext = new LocalContext(ref, vc, likelihoods, true);
        indelClassify(vc, localContext);
        isHmerIndel(vc, localContext);
        
        // Get best alleles for each read using likelihoods
        Collection<AlleleLikelihoods<GATKRead, Allele>.BestAllele> bestAlleles = likelihoods.bestAllelesBreakingTies(sampleName);
        
        // Convert to Map for easier lookup
        Map<GATKRead, Allele> readToBestAllele = new HashMap<>();
        for (AlleleLikelihoods<GATKRead, Allele>.BestAllele bestAllele : bestAlleles) {
            readToBestAllele.put(bestAllele.evidence, bestAllele.allele);
        }
        
        // Initialize bias counters for each alternate allele
        Map<Allele, int[]> biasCounts = new HashMap<>();
        for (Allele allele : vc.getAlternateAlleles()) {
            biasCounts.put(allele, new int[]{0, 0}); // [down, up]
        }
        
        // Analyze each FlowBasedPileupElement
        // Use iterator to access pileup elements
        for (PileupElement pileupElement : flowPileup) {
            FlowBasedRead flowRead = (FlowBasedRead) pileupElement.getRead();
            FlowBasedPileupElement element =
                new FlowBasedPileupElement(flowRead, pileupElement);
            
            // Get the read for this element
            GATKRead read = pileupElement.getRead();
            
            // Get best allele for this read
            Allele bestAllele = readToBestAllele.get(read);
            if (bestAllele == null || bestAllele.isReference()) {
                continue; // Skip if no best allele or reference
            }
            
            // Get call probabilities and determine bias direction
            double[] callProbs = element.getCallProbs();
            if (callProbs == null || callProbs.length == 0) {
                continue; // Skip if no probabilities available
            }
            
            int hpolCall = element.getHpolCall();
            
            // Calculate probability mass below and above the call
            double probBelow = 0.0;
            double probAbove = 0.0;
            
            for (int j = 0; j < callProbs.length; j++) {
                if (j < hpolCall) {
                    probBelow += callProbs[j];
                } else if (j > hpolCall) {
                    probAbove += callProbs[j];
                }
            }
            
            // Determine bias direction (any difference)
            boolean tendsDown = probBelow > probAbove;
            
            // Accumulate counts
            int[] counts = biasCounts.get(bestAllele);
            if (counts != null) {
                if (tendsDown) {
                    counts[0]++; // down count
                } else {
                    counts[1]++; // up count
                }
            }
        }
        
        // Format output as "|" separated tuples
        List<String> alleleBiasStrings = new ArrayList<>();
        for (Allele altAllele : vc.getAlternateAlleles()) {
            int[] counts = biasCounts.get(altAllele);
            alleleBiasStrings.add(counts[0] + "," + counts[1]);
        }
        
        String biasAnnotation = String.join("|", alleleBiasStrings);
        gb.attribute(getKeyNames().get(0), biasAnnotation);
    }
}
