package org.ultimagen.flowBasedRead.alignment;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.ultimagen.flowBasedRead.read.FlowBasedHaplotype;
import org.ultimagen.flowBasedRead.read.FlowBasedRead;
import org.ultimagen.flowBasedRead.utils.Direction;
import org.ultimagen.flowBasedRead.utils.FlowBasedAlignmentArgumentCollection;

import java.util.*;
import java.util.function.ToDoubleFunction;

/* Flow based replacement for PairHMM likelihood calculation. Likelihood calculation all-vs-all flow based reads and haplotypes

 */
public class FlowBasedAlignmentEngine implements ReadLikelihoodCalculationEngine {
    private double log10globalReadMismappingRate;
    private final double expectedErrorRatePerBase;
    private static final int ALIGNMENT_UNCERTAINTY = 4;
    final FlowBasedAlignmentArgumentCollection fbargs;
    private final Logger logger = LogManager.getLogger(this.getClass());

    /**
     * Default constructor
     * @param flowBasedArgs - arguments
     * @param log10globalReadMismappingRate - probability for wrong mapping (maximal contribution of the read to data likelihood)
     * @param expectedErrorRatePerBase - the expected rate of random sequencing errors for a read originating from its true haplotype.
     */
    public FlowBasedAlignmentEngine(final FlowBasedAlignmentArgumentCollection flowBasedArgs, final double log10globalReadMismappingRate, final double expectedErrorRatePerBase) {
        this.fbargs = flowBasedArgs;
        this.log10globalReadMismappingRate = log10globalReadMismappingRate;
        this.expectedErrorRatePerBase = expectedErrorRatePerBase;

    }

    /**
     * Read/haplotype likelihood calculation for all samples
     * @param assemblyResultSet the input assembly results.
     * @param samples the list of targeted samples.
     * @param perSampleReadList the input read sets stratified per sample.
     *
     * @param filterPoorly - if the poorly modeled reads should be removed
     * @return
     */
    @Override
    public AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods(final AssemblyResultSet assemblyResultSet,
                                                                         final SampleList samples,
                                                                         final Map<String, List<GATKRead>> perSampleReadList, final boolean filterPoorly) {
        Utils.nonNull(assemblyResultSet, "assemblyResultSet is null");
        Utils.nonNull(samples, "samples is null");
        Utils.nonNull(perSampleReadList, "perSampleReadList is null");

        final List<Haplotype> haplotypeList = assemblyResultSet.getHaplotypeList();
        final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);

        // Add likelihoods for each sample's reads to our result
        final AlleleLikelihoods<GATKRead, Haplotype> result = new AlleleLikelihoods<>(samples, haplotypes, perSampleReadList);
        final int sampleCount = result.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            final SAMFileHeader hdr= assemblyResultSet.getRegionForGenotyping().getHeader();
            computeReadLikelihoods(result.sampleMatrix(i), hdr);
        }


        result.normalizeLikelihoods(log10globalReadMismappingRate);
        if ( filterPoorly )
            result.filterPoorlyModeledEvidence(log10MinTrueLikelihood(expectedErrorRatePerBase));

        return result;
    }

    /** Calculate minimal likelihood that reasonably matching read can get. We divide errors into "expected", e.g.
     * hmer indels without 0->1/1->0 errors. Those on average occur with expectedErrorRate and "catastrophic' e.g.
     * 1->0 and 0->1 errors (or large hmer indels). Those occur with probability fbargs.fillingValue.
     * If the read has more than 3 expected and more than 2 "catastrophic" errors it will probably be deemed unfit to the
     * haplotype
     *
     * @param expectedErrorRate  error rate for expected errors.
     * @return minimal likelihood for the read to be considered not poorly modeled
     */
    private ToDoubleFunction<GATKRead> log10MinTrueLikelihood(final double expectedErrorRate) {
        final double log10ErrorRate = Math.log10(expectedErrorRate);
        final double catastrophicErrorRate = Math.log10(fbargs.fillingValue);

        return read -> {
            final double maxErrorsForRead = Math.max(3.0, Math.ceil(read.getLength() * expectedErrorRate));
            final double maxCatastrophicErrorsForRead = Math.max(2.0, Math.ceil(read.getLength() * catastrophicErrorRate));
            return maxErrorsForRead * log10ErrorRate + maxCatastrophicErrorsForRead*catastrophicErrorRate;
        };
    }


    /**
     * Compute read likelihoods for a single sample
     * @param likelihoods Single sample likelihood matrix
     * @param hdr SAM header that corresponds to the sample
     */
    private void computeReadLikelihoods(final LikelihoodMatrix<GATKRead, Haplotype> likelihoods,
                                        final SAMFileHeader hdr) {

        final List<FlowBasedRead> processedReads = new ArrayList<>(likelihoods.evidenceCount());
        final List<FlowBasedHaplotype> processedHaplotypes = new ArrayList<>(likelihoods.numberOfAlleles());
        String flowOrder = null;
        String originalFlowOrder = null;
        String fo;
        int max_class ;
        //convert all reads to FlowBasedReads (i.e. parse the matrix of P(call | haplotype) for each read from the BAM)
        for (int i = 0 ; i < likelihoods.evidenceCount(); i++) {
            final GATKRead rd = likelihoods.evidence().get(i);
            final String mc_string = hdr.getReadGroup(rd.getReadGroup()).getAttribute("mc");
            if (mc_string==null) {
                max_class = 12;
            } else {
                max_class = Integer.parseInt(mc_string);
            }

            //get flow order for conversion.
            fo = hdr.getReadGroup(rd.getReadGroup()).getFlowOrder();
            if ( fo == null ) {
                throw new GATKException("Unable to perform flow based alignment without the flow order information");
            }

            originalFlowOrder = fo.substring(0,fbargs.flowOrderCycleLength);
            final FlowBasedRead tmp = new FlowBasedRead(rd, fo, max_class, fbargs);
            tmp.applyAlignment();

            if ( flowOrder == null)  {
                fo = tmp.getFlowOrder();
                if (fo.length()>=fbargs.flowOrderCycleLength) {
                    flowOrder = fo.substring(0,fbargs.flowOrderCycleLength);
               }
            }
            processedReads.add(tmp);
        }

        //same for the haplotypes - each haplotype is converted to FlowBasedHaplotype
        FlowBasedHaplotype fbh;

        if ( flowOrder == null && hdr.getReadGroups().size() > 0 ) {
            //find the right flow order for the haplotypes (should fit to that of the reads)
            for ( final SAMReadGroupRecord rg : hdr.getReadGroups() ) {
                flowOrder = rg.getAttribute("FO");
                if ( flowOrder != null && flowOrder.length() >= fbargs.flowOrderCycleLength ) {
                    flowOrder = flowOrder.substring(0, fbargs.flowOrderCycleLength);
                    break;
                }
            }
        }

        if ( flowOrder == null ) {
            throw new GATKException("Unable to perform flow based alignment without the flow order");
        }
        for (int i = 0; i < likelihoods.numberOfAlleles(); i++){
            fbh = new FlowBasedHaplotype(likelihoods.alleles().get(i),
                        originalFlowOrder != null ? originalFlowOrder : flowOrder);
            processedHaplotypes.add(fbh);
        }

        //NOTE: we assume all haplotypes start and end on the same place!
        final int haplotypeStart = processedHaplotypes.get(0).getStart();
        final int haplotypeEnd = processedHaplotypes.get(0).getEnd();
        for (int i = 0 ; i < processedReads.size(); i++) {
            final FlowBasedRead fbr=processedReads.get(i);
            final int readStart = fbr.getStart();
            final int readEnd = fbr.getEnd();
            final int diffLeft = haplotypeStart - readStart;
            final int diffRight = readEnd - haplotypeEnd;
            //It is rare that this function is applied, maybe just some boundary cases
            //in general reads are already trimmed to the haplotype starts and ends so diff_left <= 0 and diff_right <= 0
            fbr.applyBaseClipping(Math.max(0, diffLeft), Math.max(diffRight, 0));
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


    /**
     * Aligns single read to a single haplotype. The alignment process is based on the assumption that the errors are only
     * changes in the flow calls. Thus the length of the haplotype in flow space should be equal to the length of the read
     * in the flow space (on the region where the read aligns).
     *
     * Example: align read ACGGT
     * to haplotype       GACG-TA. In flow space we see that the overlapping sequences ACGGT and ACGT are of the same length
     * and we just need to calculate P(R|H) = \prod P(Ri | Hi) for ith flow. Note that the flow matrix of the flow based
     * read {@link FlowBasedRead} allows for trivial calculation of this product.
     *
     * As an additional optimization we note that almost always the interval to be aligned corresponds to the overlap
     * of the read and the haplotype to the reference. We allow for some uncertainty in this, especially when the read
     * falls inside the deletion in the haplotype relative to the reference.
     *
     * @param haplotype FlowBasedHaplotype single haplotype
     * @param read FlowBasedRead single read (trimmed to the haplotype)
     * @return
     * @throws GATKException
     */
    private double haplotypeReadMatching(final FlowBasedHaplotype haplotype, final FlowBasedRead read) throws GATKException {

        if (read.getDirection() != Direction.REFERENCE ) {
            throw new GATKException.ShouldNeverReachHereException("Read should be aligned with the reference");
        }

        if (!read.isTrimmedToHaplotype()) {
            throw new GATKException.ShouldNeverReachHereException("Reads should be trimmed to the haplotype");
        }

        if (!read.isValid()) {
            return Double.NEGATIVE_INFINITY;
        }

        // the read is assumed to be trimmed to the haplotype by ReadThreadingAssembler.finalizeRegion, the region of the
        // haplotype to be aligned is estimated by finding the points on the haplotypes that align to the start and the end
        // of the read on  the reference
        final int haplotypeStart = ReadUtils.getReadIndexForReferenceCoordinate(haplotype.getStart(), haplotype.getCigar(),
                read.getTrimmedStart()).getLeft();

        final int haplotypeEnd = ReadUtils.getReadIndexForReferenceCoordinate(haplotype.getStart(), haplotype.getCigar(),
                read.getTrimmedEnd()).getLeft();

        final int haplotypeLength = haplotypeEnd - haplotypeStart;
        final int readLength = read.seqLength();


        //in case there is a deletion on the haplotype and hte read falls inside the deletion (thus length of the read is
        // longer than the length of the trimmed haplotype:
        // Reference:  GACACACACACT
        // Read:              ACACT
        // Haplotype   GACAC------T
        // In this case the read (that is not informative) is aligned inside the deletion and thus if we clip the haplotype
        // to the start position of the read will not support it (haplotype: T, read: ACACT), so in cases we see that the resulting
        // haploype length is shorter than the resulting read length we extend the "trimmed haplotype" by this length and
        // also allow more uncertainty for the actual starting position of the alignment
        final int uncertainty = Math.max(readLength - haplotypeLength,0);
        final int leftClip = Math.max(haplotypeStart-uncertainty,0);
        final int rightClip = Math.max(haplotype.length()-haplotypeEnd-1-uncertainty,0);


        if ((leftClip < 0) || (rightClip < 0)  || (leftClip >= haplotype.length() ) || ( rightClip >= haplotype.length())) {
            return 1;
        }

        int [] tmp = haplotype.findLeftClipping(leftClip);
        int clipLeft = tmp[0];
        final int leftHmerClip = tmp[1];

        tmp = haplotype.findRightClipping(rightClip);
        int clipRight = tmp[0];
        final int rightHmerClip = tmp[1];


        if ((clipLeft >= haplotype.getKeyLength()) || (clipRight >= haplotype.getKeyLength())){
            return Double.NEGATIVE_INFINITY;
        }

        if ((leftHmerClip <0) | (rightHmerClip < 0)) {
            throw new GATKException.ShouldNeverReachHereException("Negative hmer clips found. Check");
        }

        final int originalLength = haplotype.getKeyLength();
        clipLeft = clipLeft-ALIGNMENT_UNCERTAINTY>0 ? clipLeft-ALIGNMENT_UNCERTAINTY:0;
        clipRight = originalLength - clipRight + ALIGNMENT_UNCERTAINTY < originalLength ?
                            originalLength - clipRight+ALIGNMENT_UNCERTAINTY : originalLength;

        final int [] key;
        key = Arrays.copyOfRange(haplotype.getKey(), clipLeft, clipRight);

        final byte[] flowOrder;
        flowOrder = Arrays.copyOfRange(haplotype.getFlowOrderArray(), clipLeft, clipRight);
        int startingPoint = 0;
        for (int i = 0; i < flowOrder.length; i++) {
            if (flowOrder[i] == read.getFlowOrderArray()[0]) {
                startingPoint = i;
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
            logger.debug("starting point: " + read.getName() + " " + startingPoint);
        }

        // this is the heart of the calculation of the likelihood. Since we assume in our model that between the
        // read and the haplotype there are no insertions / deletions of flows, we just search for the best starting
        // position of the read on the haplotype and then to calculate P(R|H) we just select the corresponding probabilities
        // from the flow matrix.
        // Another important optimization is that the best starting position of the read on the haplotype can be calculated
        // from the starting position of the read on the reference and the alignment of the haplotype to the reference.
        // This is why we trim the haplotype so that the aligned intervals have the same length as the read and allow a
        // small uncertainty.
        double bestAlignment = Double.NEGATIVE_INFINITY;
        for ( int s = startingPoint ; ( s+read.getKeyLength() <= key.length); s+= ALIGNMENT_UNCERTAINTY){
            final int [] locationsToFetch = new int[read.getKeyLength()];
            for (int i = s; i < s+read.getKeyLength(); i++ ){
                locationsToFetch[i-s] = key[i]&0xff;
                locationsToFetch[i-s] = locationsToFetch[i-s] < read.getMaxHmer()+1 ?
                        locationsToFetch[i-s] : read.getMaxHmer()+1;
            }
            double result = 0 ;
            for ( int i = 0 ; i < locationsToFetch.length; i++ ){
                final double prob;
                result += Math.log10(prob = read.getProb(i, locationsToFetch[i]));
                if ( logger.isDebugEnabled() )
                    logger.debug("prob:" + read.getName() + " " + i + " " + locationsToFetch[i] + " " + prob);
            }
            if (result > bestAlignment) {
                bestAlignment = result;
            }
        }
        return bestAlignment;
    }

    //This function is for testing purposes only
    @VisibleForTesting
    public AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods(
            final List<Haplotype> haplotypeList, final List<GATKRead> reads, final boolean filterPoorly) {


        final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        final ArrayList<String> _sampList = new ArrayList<>();
        _sampList.add("HG001");
        final SampleList samples = new IndexedSampleList(_sampList);

        // Add likelihoods for each sample's reads to our result
        final HashMap<String, List<GATKRead>> perSampleReadList = new HashMap<>();
        perSampleReadList.put("HG001", reads);

        final AlleleLikelihoods<GATKRead, Haplotype> result = new AlleleLikelihoods<>(samples, haplotypes, perSampleReadList);
        final int sampleCount = result.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            computeReadLikelihoods(result.sampleMatrix(i), null);
        }

        result.normalizeLikelihoods(log10globalReadMismappingRate);
        if ( filterPoorly )
            result.filterPoorlyModeledEvidence(log10MinTrueLikelihood(expectedErrorRatePerBase));

        return result;
    }

    @Override
    public void close() {}

}

