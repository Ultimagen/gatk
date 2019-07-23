package org.ultimagenomics.flow_based_read.alignment;

import org.ultimagenomics.flow_based_read.utils.Direction;
import org.ultimagenomics.flow_based_read.read.FlowBasedHaplotype;
import org.ultimagenomics.flow_based_read.read.FlowBasedRead;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.*;


public class FlowBasedAlignmentEngine implements ReadLikelihoodCalculationEngine {
    final private static String FLOW_ORDER="TACG";
    private double log10globalReadMismappingRate;
    private static final double EXPECTED_ERROR_RATE_PER_BASE = 0.02;

    public FlowBasedAlignmentEngine(double log10globalReadMismappingRate) {
        this.log10globalReadMismappingRate = log10globalReadMismappingRate;

    }
    @Override
    public ReadLikelihoods<Haplotype> computeReadLikelihoods(AssemblyResultSet assemblyResultSet, SampleList samples, Map<String, List<GATKRead>> perSampleReadList) {
        Utils.nonNull(assemblyResultSet, "assemblyResultSet is null");
        Utils.nonNull(samples, "samples is null");
        Utils.nonNull(perSampleReadList, "perSampleReadList is null");

        final List<Haplotype> haplotypeList = assemblyResultSet.getHaplotypeList();
        final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);

        // Add likelihoods for each sample's reads to our result
        final ReadLikelihoods<Haplotype> result = new ReadLikelihoods<>(samples, haplotypes, perSampleReadList);
        final int sampleCount = result.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            computeReadLikelihoods(result.sampleMatrix(i));
        }

        result.normalizeLikelihoods(log10globalReadMismappingRate);
        result.filterPoorlyModeledReads(EXPECTED_ERROR_RATE_PER_BASE);

        return result;
    }

    public ReadLikelihoods<Haplotype> computeReadLikelihoods(final List<Haplotype> haplotypeList,
                                                             final List<GATKRead> reads) {


        final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        ArrayList<String> _sampList = new ArrayList<>();
        _sampList.add("HG001");
        SampleList samples = new IndexedSampleList(_sampList);

        // Add likelihoods for each sample's reads to our result
        HashMap<String, List<GATKRead>> perSampleReadList = new HashMap<>();
        perSampleReadList.put("HG001", reads);

        final ReadLikelihoods<Haplotype> result = new ReadLikelihoods<>(samples, haplotypes, perSampleReadList);
        final int sampleCount = result.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            computeReadLikelihoods(result.sampleMatrix(i));
        }

        result.normalizeLikelihoods(log10globalReadMismappingRate);
        result.filterPoorlyModeledReads(EXPECTED_ERROR_RATE_PER_BASE);

        return result;
    }


    private void computeReadLikelihoods(LikelihoodMatrix<Haplotype> likelihoods) {

        List<FlowBasedRead> processedReads = new ArrayList<>(likelihoods.numberOfReads());
        for (int i = 0 ; i < likelihoods.numberOfReads(); i++) {
            FlowBasedRead tmp = new FlowBasedRead(likelihoods.reads().get(i));
            tmp.apply_alignment();
            processedReads.add(tmp);
        }

        FlowBasedHaplotype fbh;

        for (int i = 0; i < likelihoods.numberOfAlleles(); i++){
            fbh = new FlowBasedHaplotype(likelihoods.alleles().get(i), FLOW_ORDER, 8);
            for (int j = 0 ; j < likelihoods.numberOfReads(); j++){

                likelihoods.set(i,j,haplotypeReadMatching(fbh,processedReads.get(j)));
            }
        }
    }

    private double haplotypeReadMatching(FlowBasedHaplotype haplotype, FlowBasedRead read) throws GATKException {
        if (read.getDirection() != Direction.REFERENCE ) {
            throw new GATKException.ShouldNeverReachHereException("Read should be aligned with the reference");
        }
        if ((haplotype.getStart()>read.getStart())||(haplotype.getEnd() < read.getEnd())) {
            throw new GATKException.ShouldNeverReachHereException("Read should be contained in the haplotype");
        }

        int read_start = read.getStart();
        int read_end = read.getEnd();
        int haplotype_start = ReadUtils.getReadCoordinateForReferenceCoordinate(haplotype.getStart(), haplotype.getCigar(), read_start, null, false);
        int haplotype_end = ReadUtils.getReadCoordinateForReferenceCoordinate(haplotype.getStart(), haplotype.getCigar(), read_end, null, false);

        int left_clip = haplotype_start;
        int right_clip = haplotype.length()-haplotype_end-1;


        if ((left_clip < 0) || (right_clip < 0)  || (left_clip >= haplotype.length() ) || ( right_clip >= haplotype.length())) {
            return 1;
        }

        int [] tmp = haplotype.find_left_clipping(left_clip);
        int clip_left = tmp[0];
        int left_hmer_clip = tmp[1];

        tmp = haplotype.find_right_clipping(right_clip);
        int clip_right = tmp[0];
        int right_hmer_clip = tmp[1];

        if ((left_hmer_clip > 11) || (right_hmer_clip > 11)) {
            throw new GATKException.ShouldNeverReachHereException("Weird haplotype clip calculated");
        }

        if ((clip_left >= haplotype.getKeyLength()) || (clip_right >= haplotype.getKeyLength())){
            return Double.NEGATIVE_INFINITY;
        }


        if ((left_hmer_clip <0) | (right_hmer_clip < 0)) {
            throw new GATKException.ShouldNeverReachHereException("Negative hmer clips found. Check");
        }

        int original_length = haplotype.getKeyLength();
        clip_left = clip_left-4>0 ? clip_left-4:0;
        clip_right = original_length - clip_right + 4 < original_length ? original_length - clip_right+4 : original_length;

        byte [] key;
        key = Arrays.copyOfRange(haplotype.getKey(), clip_left, clip_right);

        byte[] flow_order;
        flow_order = Arrays.copyOfRange(haplotype.getFlowOrderArray(), clip_left, clip_right);
        int starting_point = 0;
        for (int i = 0; i < flow_order.length; i++) {
            if (flow_order[i] == read.getFlowOrderArray()[0]) {
                starting_point = i;
                break;
            }
        }
        double best_alignment = Double.NEGATIVE_INFINITY;
        for ( int s = starting_point ; (s < starting_point + FLOW_ORDER.length()*2) && ( s+read.getKeyLength() <= key.length); s+= FLOW_ORDER.length()){
            int [] locations_to_fetch = new int[read.getKeyLength()];
            for (int i = s; i < s+read.getKeyLength(); i++ ){
                locations_to_fetch[i-s] = key[i]&0xff;
            }
            double result = 0 ;
            for ( int i = 0 ; i < locations_to_fetch.length; i++ ){
                result += Math.log10(read.getProb(i, locations_to_fetch[i]));
            }
            if (result > best_alignment) {
                best_alignment = result;
            }
        }
        return best_alignment;
    }


    @Override
    public void close() {}

}

