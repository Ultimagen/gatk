package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.FlowBasedAlignmentArgumentCollection;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.*;
import org.broadinstitute.hellbender.utils.haplotype.FlowBasedHaplotype;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.function.ToDoubleFunction;
import java.util.stream.IntStream;

/* Flow based replacement for PairHMM likelihood calculation. Likelihood calculation all-vs-all flow based reads and haplotypes

 */
public class FlowBasedAlignmentEngine implements ReadLikelihoodCalculationEngine {
    public static final double MAX_ERRORS_FOR_READ_CAP = 3.0;
    public static final double MAX_CATASTROPHIC_ERRORS_FOR_READ_CAP = 2.0;

    private double log10globalReadMismappingRate;
    private final double expectedErrorRatePerBase;
    private final boolean symmetricallyNormalizeAllelesToReference;

    private static final int ALIGNMENT_UNCERTAINTY = 4;
    final FlowBasedAlignmentArgumentCollection fbargs;
    private final Logger logger = LogManager.getLogger(this.getClass());

    private final double commonProbValue = 0.001;
    private final boolean dynamicReadDisqualification;
    private final double readDisqualificationScale;

    private ForkJoinPool    threadPool;
    static final double prob0 = 0.988;
    static final double prob0log10 = Math.log10(prob0);
    static final double prob1 = 0.001;
    static final double prob1log10 = Math.log10(prob1);


    /**
     * Default constructor
     * @param fbargs - arguments
     * @param log10globalReadMismappingRate - probability for wrong mapping (maximal contribution of the read to data likelihood)
     * @param expectedErrorRatePerBase - the expected rate of random sequencing errors for a read originating from its true haplotype.
     */
    public FlowBasedAlignmentEngine(final FlowBasedAlignmentArgumentCollection fbargs, final double log10globalReadMismappingRate, final double expectedErrorRatePerBase, final boolean dynamicReadDisqualification, final double readDisqualificationScale) {
        this.fbargs = fbargs;
        this.log10globalReadMismappingRate = log10globalReadMismappingRate;
        this.expectedErrorRatePerBase = expectedErrorRatePerBase;
        this.readDisqualificationScale = readDisqualificationScale;
        this.dynamicReadDisqualification = dynamicReadDisqualification;
        this.symmetricallyNormalizeAllelesToReference = true;

        if ( this.fbargs.flowLikelihoodParallelThreads > 0 ) {
            threadPool = new ForkJoinPool(this.fbargs.flowLikelihoodParallelThreads);
        }
    }

    /**
     * Read/haplotype likelihood calculation for all samples
     * @param samples the list of targeted samples.
     * @param perSampleReadList the input read sets stratified per sample.
     *
     * @param filterPoorly - if the poorly modeled reads should be removed
     * @return
     */
    public AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods(final List<Haplotype> haplotypeList,
                                                                         final SAMFileHeader hdr,
                                                                         final SampleList samples,
                                                                         final Map<String, List<GATKRead>> perSampleReadList, final boolean filterPoorly) {
        Utils.nonNull(samples, "samples is null");
        Utils.nonNull(perSampleReadList, "perSampleReadList is null");
        Utils.nonNull(haplotypeList, "haplotypeList is null");

        final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);

        // Add likelihoods for each sample's reads to our result
        final AlleleLikelihoods<GATKRead, Haplotype> result = new AlleleLikelihoods<>(samples, haplotypes, perSampleReadList);
        final int sampleCount = result.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            computeReadLikelihoods(result.sampleMatrix(i), hdr);
        }

        result.normalizeLikelihoods(log10globalReadMismappingRate, symmetricallyNormalizeAllelesToReference);
        filterPoorlyModeledEvidence(result, dynamicReadDisqualification, expectedErrorRatePerBase, readDisqualificationScale);

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
    @Override
    public ToDoubleFunction<GATKRead> log10MinTrueLikelihood(final double expectedErrorRate, final boolean capLikelihoods) {
        final double log10ErrorRate = Math.log10(expectedErrorRate);
        final double catastrophicErrorRate = Math.log10(fbargs.fillingValue);

        return read -> {
            final double maxErrorsForRead = capLikelihoods ? Math.max(MAX_ERRORS_FOR_READ_CAP, Math.ceil(read.getLength() * expectedErrorRate)) : Math.ceil(read.getLength() * expectedErrorRate);
            final double maxCatastrophicErrorsForRead = capLikelihoods ? Math.max(MAX_CATASTROPHIC_ERRORS_FOR_READ_CAP, Math.ceil(read.getLength() * catastrophicErrorRate)) : Math.ceil(read.getLength() * catastrophicErrorRate);
            return maxErrorsForRead * log10ErrorRate + maxCatastrophicErrorsForRead * catastrophicErrorRate;
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

        // establish flow order based on the first evidence. Note that all reads belong to the same sample (group)
        final FlowBasedReadUtils.ReadGroupInfo rgInfo = (likelihoods.evidenceCount() != 0)
                ? FlowBasedReadUtils.getReadGroupInfo(hdr, likelihoods.evidence().get(0))
                : null;
        final String flowOrder = (rgInfo != null)
                ? rgInfo.flowOrder.substring(0, fbargs.flowOrderCycleLength)
                : FlowBasedReadUtils.findFirstUsableFlowOrder(hdr, fbargs);

        //convert all reads to FlowBasedReads (i.e. parse the matrix of P(call | haplotype) for each read from the BAM)
        for (int i = 0 ; i < likelihoods.evidenceCount(); i++) {
            final GATKRead rd = likelihoods.evidence().get(i);

            final FlowBasedRead fbRead = new FlowBasedRead(rd, flowOrder, rgInfo.maxClass, fbargs);
            fbRead.applyAlignment();

            processedReads.add(fbRead);
        }

        if ( flowOrder == null ) {
            throw new GATKException("Unable to perform flow based alignment without the flow order");
        }

        for (int i = 0; i < likelihoods.numberOfAlleles(); i++){
            final FlowBasedHaplotype fbh = new FlowBasedHaplotype(likelihoods.alleles().get(i), flowOrder);
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
            fbr.applyBaseClipping(Math.max(0, diffLeft), Math.max(diffRight, 0), true);
        }

        for (int i = 0; i < likelihoods.numberOfAlleles(); i++){
            final FlowBasedHaplotype  fbh = processedHaplotypes.get(i);
            if ( threadPool != null  ) {
                final int haplotypeIndex = i;
                Callable<Void>  callable = () -> {
                    IntStream.range(0, likelihoods.evidenceCount()).parallel().forEach(j -> {
                        final double likelihood = haplotypeReadMatching(fbh, processedReads.get(j));
                        likelihoods.set(haplotypeIndex, j, likelihood);
                    });
                    return null;
                };
                try {
                    threadPool.submit(callable).get();
                } catch (InterruptedException | ExecutionException e) {
                    throw new RuntimeException(e);
                }
            } else {
                for (int j = 0; j < likelihoods.evidenceCount(); j++) {
                    final double likelihood = haplotypeReadMatching(fbh, processedReads.get(j));
                    likelihoods.set(i, j, likelihood);
                }
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
    public double haplotypeReadMatching(final FlowBasedHaplotype haplotype, final FlowBasedRead read) throws GATKException {

        if (read.getDirection() != FlowBasedRead.Direction.REFERENCE ) {
            throw new GATKException.ShouldNeverReachHereException("Read should be aligned with the reference");
        }

        if (!read.isBaseClipped()) {
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

        int [] leftClipping = haplotype.findLeftClipping(leftClip);
        int clipLeft = leftClipping[0];
        final int leftHmerClip = leftClipping[1];

        leftClipping = haplotype.findRightClipping(rightClip);
        int clipRight = leftClipping[0];
        final int rightHmerClip = leftClipping[1];


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

        if ( fbargs.flowLikelihoodOptimizedComp ) {

            int[]   key_o = haplotype.getKey();
            int     key_ostart = clipLeft;
            int     key_olen = clipRight - clipLeft;

            byte[]  flow_order_o = haplotype.getFlowOrderArray();
            int     flow_order_ostart = clipLeft;
            int     flow_order_olen = clipRight - clipLeft;
            byte    read_flow_0 = read.getFlowOrderArray()[0];
            int startingPoint = 0;
            for (int i = 0; i < flow_order_olen; i++) {
                if (flow_order_o[i+flow_order_ostart] == read_flow_0) {
                    startingPoint = i;
                    break;
                }
            }

            if (logger.isDebugEnabled()) {
                logger.debug("haplotype   " + read.getName() + " " + haplotype.getBaseString());
                logger.debug("read        " + read.getName() + " " + read.getBasesString());
                logger.debug("haplotype.f " + read.getName() + " " + new String(haplotype.getFlowOrderArray()));
                logger.debug("read.f      " + read.getName() + " " + new String(read.getFlowOrderArray()));
                logger.debug("haplotype.K " + read.getName() + " " + FlowBasedKeyCodec.keyAsString(haplotype.getKey()));
                logger.debug("haplotype.k " + read.getName() + " " + FlowBasedKeyCodec.keyAsString(Arrays.copyOfRange(haplotype.getKey(), clipLeft, clipRight)));
                logger.debug("read.k      " + read.getName() + " " + FlowBasedKeyCodec.keyAsString(read.getKey()));
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
            int    read_maxhmer_1 = read.getMaxHmer() + 1;
            int    read_key_len = read.getKeyLength();
            final int[] locationsToFetch = new int[read_key_len];
            for (int s = startingPoint; (s + read_key_len <= key_olen); s += ALIGNMENT_UNCERTAINTY) {
                for (int i = s; i < s + read_key_len; i++) {
                    if ( (locationsToFetch[i - s] = key_o[i+key_ostart] & 0xff) >= read_maxhmer_1 )
                        locationsToFetch[i - s] = read_maxhmer_1;
                }
                double result = 0;
                for (int i = 0; i < read_key_len; i++) {
                    final double prob = read.getProb(i, locationsToFetch[i]);
                    double      log10;
                    if ( prob == prob0 )
                        log10 = prob0log10;
                    else if ( prob == prob1 )
                        log10 = prob1log10;
                    else
                        log10 = Math.log10(prob);
                    result += log10;
                    if (logger.isDebugEnabled())
                        logger.debug("prob:" + read.getName() + " " + i + " " + locationsToFetch[i] + " " + prob);

                    // shortcut - result will never rise
                    if ( result < bestAlignment )
                        break;
                }
                if (result > bestAlignment) {
                    bestAlignment = result;
                }
            }
            return bestAlignment;

        } else {

            final int[] key;
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

            if (logger.isDebugEnabled()) {
                logger.debug("haplotype   " + read.getName() + " " + haplotype.getBaseString());
                logger.debug("read        " + read.getName() + " " + read.getBasesString());
                logger.debug("haplotype.f " + read.getName() + " " + new String(haplotype.getFlowOrderArray()));
                logger.debug("read.f      " + read.getName() + " " + new String(read.getFlowOrderArray()));
                logger.debug("haplotype.K " + read.getName() + " " + FlowBasedKeyCodec.keyAsString(haplotype.getKey()));
                logger.debug("haplotype.k " + read.getName() + " " + FlowBasedKeyCodec.keyAsString(key));
                logger.debug("read.k      " + read.getName() + " " + FlowBasedKeyCodec.keyAsString(read.getKey()));
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
            for (int s = startingPoint; (s + read.getKeyLength() <= key.length); s += ALIGNMENT_UNCERTAINTY) {
                final int[] locationsToFetch = new int[read.getKeyLength()];
                for (int i = s; i < s + read.getKeyLength(); i++) {
                    locationsToFetch[i - s] = key[i] & 0xff;
                    locationsToFetch[i - s] = locationsToFetch[i - s] < read.getMaxHmer() + 1 ?
                            locationsToFetch[i - s] : read.getMaxHmer() + 1;
                }
                double result = 0;
                for (int i = 0; i < locationsToFetch.length; i++) {
                    final double prob;
                    result += Math.log10(prob = read.getProb(i, locationsToFetch[i]));
                    if (logger.isDebugEnabled()) {
                        logger.debug("prob:" + read.getName() + " " + i + " " + locationsToFetch[i] + " " + prob);
                    }
                }
                if (result > bestAlignment) {
                    bestAlignment = result;
                }
            }
            return bestAlignment;
        }
    }

    //This function is for testing purposes only
    @VisibleForTesting
    public AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods(
            final List<Haplotype> haplotypeList, final List<GATKRead> reads, final boolean filterPoorly, final SAMFileHeader hdr) {


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
            computeReadLikelihoods(result.sampleMatrix(i), hdr);
        }

        result.normalizeLikelihoods(log10globalReadMismappingRate, symmetricallyNormalizeAllelesToReference);
        if ( filterPoorly ) {
            result.filterPoorlyModeledEvidence(log10MinTrueLikelihood(expectedErrorRatePerBase, false));
        }

        return result;
    }

    @Override
    public void close() {}

}

