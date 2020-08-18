package org.ultimagenomics.flow_based_read.alignment;

import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
import org.ultimagenomics.flow_based_read.utils.FlowBasedAlignmentArgumentCollection;

import java.util.*;
import java.util.function.ToDoubleFunction;


public class FlowBasedAlignmentEngine implements ReadLikelihoodCalculationEngine {
    private double log10globalReadMismappingRate;
    private final double expectedErrorRatePerBase;
    private static final int ALIGNMENT_UNCERTAINTY = 4;
    private static final double LOG10_QUAL_PER_BASE = Double.NEGATIVE_INFINITY;
    final FlowBasedAlignmentArgumentCollection fbargs;
    private final Logger logger = LogManager.getLogger(this.getClass());

    public FlowBasedAlignmentEngine(final FlowBasedAlignmentArgumentCollection flowBasedArgs, double log10globalReadMismappingRate, final double expectedErrorRatePerBase) {
        this.fbargs = flowBasedArgs;

        this.log10globalReadMismappingRate = log10globalReadMismappingRate;
        this.expectedErrorRatePerBase = expectedErrorRatePerBase;

    }
    @Override
    public AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods(AssemblyResultSet assemblyResultSet,
                                                             SampleList samples,
                                                             Map<String, List<GATKRead>> perSampleReadList) {
        Utils.nonNull(assemblyResultSet, "assemblyResultSet is null");
        Utils.nonNull(samples, "samples is null");
        Utils.nonNull(perSampleReadList, "perSampleReadList is null");

        final List<Haplotype> haplotypeList = assemblyResultSet.getHaplotypeList();
        final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);

        // Add likelihoods for each sample's reads to our result
        final AlleleLikelihoods<GATKRead, Haplotype> result = new AlleleLikelihoods<>(samples, haplotypes, perSampleReadList);
        final int sampleCount = result.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            SAMFileHeader hdr= assemblyResultSet.getRegionForGenotyping().getHeader();
            computeReadLikelihoods(result.sampleMatrix(i), hdr);
        }


        result.normalizeLikelihoods(log10globalReadMismappingRate);
        result.filterPoorlyModeledEvidence(log10MinTrueLikelihood(expectedErrorRatePerBase));

        return result;
    }

    private ToDoubleFunction<GATKRead> log10MinTrueLikelihood(final double expectedErrorRate) {
        final double log10ErrorRate = Math.log10(expectedErrorRate);
        final double catastrophicErrorRate = Math.log10(fbargs.filling_value);

        return read -> {
            final double maxErrorsForRead = Math.max(3.0, Math.ceil(read.getLength() * expectedErrorRate));
            final double maxCatastrophicErrorsForRead = Math.max(2.0, Math.ceil(read.getLength() * catastrophicErrorRate));
            return maxErrorsForRead * log10ErrorRate + maxCatastrophicErrorsForRead*catastrophicErrorRate;
        };
    }



    private void computeReadLikelihoods(LikelihoodMatrix<GATKRead, Haplotype> likelihoods,
                                        SAMFileHeader hdr) {

        List<FlowBasedRead> processedReads = new ArrayList<>(likelihoods.evidenceCount());
        List<FlowBasedHaplotype> processedHaplotypes = new ArrayList<>(likelihoods.numberOfAlleles());
        String flow_order = null;
        String original_flow_order = null;
        String fo;
        int max_class = 12 ;
        for (int i = 0 ; i < likelihoods.evidenceCount(); i++) {
            GATKRead rd = likelihoods.evidence().get(i);
            String mc_string = hdr.getReadGroup(rd.getReadGroup()).getAttribute("mc");
            if (mc_string==null) {
                max_class = 12;
            } else {
                max_class = Integer.parseInt(mc_string);
            }

            fo = hdr.getReadGroup(rd.getReadGroup()).getFlowOrder();
            original_flow_order = fo.substring(0,fbargs.flowOrderCycleLength);
            FlowBasedRead tmp = new FlowBasedRead(rd, fo, max_class, fbargs);
            tmp.apply_alignment();

            if ( flow_order == null)  {
                fo = tmp.getFlowOrder();
                if (fo.length()>=fbargs.flowOrderCycleLength) {
                    flow_order = fo.substring(0,fbargs.flowOrderCycleLength);
               }
            }
            processedReads.add(tmp);
        }

        FlowBasedHaplotype fbh;

        for (int i = 0; i < likelihoods.numberOfAlleles(); i++){
            fbh = new FlowBasedHaplotype(likelihoods.alleles().get(i), original_flow_order, max_class);
            processedHaplotypes.add(fbh);
        }

        //NOTE: we assume all haplotypes start and end on the same place!
        int haplotype_start = processedHaplotypes.get(0).getStart();
        int haplotype_end = processedHaplotypes.get(0).getEnd();
        for (int i = 0 ; i < processedReads.size(); i++) {
            FlowBasedRead fbr=processedReads.get(i);
            int read_start = fbr.getStart();
            int read_end = fbr.getEnd();
            int diff_left = haplotype_start - read_start;
            int diff_right = read_end - haplotype_end;
            fbr.apply_base_clipping(Math.max(0, diff_left), Math.max(diff_right, 0));
        }

        for (int i = 0; i < likelihoods.numberOfAlleles(); i++){
            fbh = processedHaplotypes.get(i);
            for (int j = 0 ; j < likelihoods.evidenceCount(); j++){
                final double    likelihood = haplotypeReadMatching(fbh,processedReads.get(j));
                likelihoods.set(i,j,likelihood);
                if ( logger.isDebugEnabled() )
                    logger.debug("likelihood: " + likelihood + " "  + processedReads.get(j).getName() + " " + fbh.getBaseString());
            }
        }

    }

    //This function is for testing purposes only
    public AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods(
            final List<Haplotype> haplotypeList, final List<GATKRead> reads) {


        final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        ArrayList<String> _sampList = new ArrayList<>();
        _sampList.add("HG001");
        SampleList samples = new IndexedSampleList(_sampList);

        // Add likelihoods for each sample's reads to our result
        HashMap<String, List<GATKRead>> perSampleReadList = new HashMap<>();
        perSampleReadList.put("HG001", reads);

        final AlleleLikelihoods<GATKRead, Haplotype> result = new AlleleLikelihoods<>(samples, haplotypes, perSampleReadList);
        final int sampleCount = result.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            computeReadLikelihoods(result.sampleMatrix(i), null);
        }

        result.normalizeLikelihoods(log10globalReadMismappingRate);
        result.filterPoorlyModeledEvidence(log10MinTrueLikelihood(expectedErrorRatePerBase));

        return result;
    }

    private double haplotypeReadMatching(FlowBasedHaplotype haplotype, FlowBasedRead read) throws GATKException {

        if (read.getDirection() != Direction.REFERENCE ) {
            throw new GATKException.ShouldNeverReachHereException("Read should be aligned with the reference");
        }

        if (!read.isTrimmed_to_haplotype()) {
            throw new GATKException.ShouldNeverReachHereException("Reads should be trimmed to the haplotype");
        }

        if (!read.is_valid()) return Double.NEGATIVE_INFINITY;

        int haplotype_start = ReadUtils.getReadIndexForReferenceCoordinate(haplotype.getStart(), haplotype.getCigar(),
                read.getTrimmedStart()).getLeft();
        int haplotype_end = ReadUtils.getReadIndexForReferenceCoordinate(haplotype.getStart(), haplotype.getCigar(),
                read.getTrimmedEnd()).getLeft();

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


        if ((clip_left >= haplotype.getKeyLength()) || (clip_right >= haplotype.getKeyLength())){
            return Double.NEGATIVE_INFINITY;
        }

        if ((left_hmer_clip <0) | (right_hmer_clip < 0)) {
            throw new GATKException.ShouldNeverReachHereException("Negative hmer clips found. Check");
        }

        int original_length = haplotype.getKeyLength();
        clip_left = clip_left-ALIGNMENT_UNCERTAINTY>0 ? clip_left-ALIGNMENT_UNCERTAINTY:0;
        clip_right = original_length - clip_right + ALIGNMENT_UNCERTAINTY < original_length ?
                            original_length - clip_right+ALIGNMENT_UNCERTAINTY : original_length;

        int [] key;
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

        if ( logger.isDebugEnabled() ) {
            logger.debug("haplotype   " + read.getName() + " " + haplotype.getBaseString());
            logger.debug("read        " + read.getName() + " " + read.getBasesString());
            logger.debug("haplotype.f " + read.getName() + " " + new String(haplotype.getFlowOrderArray()));
            logger.debug("read.f      " + read.getName() + " " + new String(read.getFlowOrderArray()));
            logger.debug("haplotype.K " + read.getName() + " " + FlowBasedRead.keyAsString(haplotype.getKey()));
            logger.debug("haplotype.k " + read.getName() + " " + FlowBasedRead.keyAsString(key));
            logger.debug("read.k      " + read.getName() + " " + FlowBasedRead.keyAsString(read.getKey()));
            logger.debug("starting point: " + read.getName() + " " + starting_point);
        }

        /*
        if ( logger.isDebugEnabled() ) {
            read.logMatrix(logger, "haplotypeReadMatching");
        }
         */

        double best_alignment = Double.NEGATIVE_INFINITY;
        for ( int s = starting_point ; (s < starting_point + ALIGNMENT_UNCERTAINTY*2+1) &&
                ( s+read.getKeyLength() <= key.length); s+= ALIGNMENT_UNCERTAINTY){
            int [] locations_to_fetch = new int[read.getKeyLength()];
            for (int i = s; i < s+read.getKeyLength(); i++ ){
                locations_to_fetch[i-s] = key[i]&0xff;
                locations_to_fetch[i-s] = locations_to_fetch[i-s] < read.getMaxHmer()+1 ?
                        locations_to_fetch[i-s] : read.getMaxHmer()+1;
            }
            double result = 0 ;
            for ( int i = 0 ; i < locations_to_fetch.length; i++ ){
                double prob;
                result += Math.log10(prob = read.getProb(i, locations_to_fetch[i]));
                if ( logger.isDebugEnabled() )
                    logger.debug("prob:" + read.getName() + " " + i + " " + locations_to_fetch[i] + " " + prob);
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

