package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pairhmm.FlowBasedPairHMM;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputation;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputator;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.haplotype.FlowBasedHaplotype;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.tools.FlowBasedAlignmentArgumentCollection;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.ToDoubleFunction;

/**
 *  Flow Based HMM, intended to incorporate the scoring model of the {@link FlowBasedAlignmentEngine} while allowing for frame-shift insertions
 *  and deletions for better genotyping.
 */
public class FlowBasedHMMEngine implements ReadLikelihoodCalculationEngine {
    private final double readDisqualificationScale;
    private final boolean dynamicReadDisqualification;
    private double log10globalReadMismappingRate;
    private final double expectedErrorRatePerBase;
    private PairHMMLikelihoodCalculationEngine.PCRErrorModel pcrErrorModel;
    final FlowBasedAlignmentArgumentCollection fbargs;

    @VisibleForTesting
    public static final double INITIAL_QSCORE = 40.0;

    private final FlowBasedPairHMM flowPairHMM;

    public static final byte MIN_USABLE_Q_SCORE_DEFAULT = 6;
    private static final int MIN_ADJUSTED_QSCORE = 10;

    private byte minUsableIndelScoreToUse;

    private PairHMMInputScoreImputator inputScoreImputator;
    private final DragstrParams dragstrParams;
    private byte constantGCP;
    private final byte flatInsertionPenalty;
    private final byte flatDeletionPenalty;


    /**
     * Default constructor
     * @param flowBasedArgs - arguments
     * @param log10globalReadMismappingRate - probability for wrong mapping (maximal contribution of the read to data likelihood)
     * @param expectedErrorRatePerBase - the expected rate of random sequencing errors for a read originating from its true haplotype.
     * @param pcrErrorModel
     */
    public FlowBasedHMMEngine(final FlowBasedAlignmentArgumentCollection flowBasedArgs, final byte constantGCP,  final double log10globalReadMismappingRate, final double expectedErrorRatePerBase,
                              final PairHMMLikelihoodCalculationEngine.PCRErrorModel pcrErrorModel, final DragstrParams dragstrParams, final boolean dynamicReadDisqualification, final double readDisqualificationScale,
                              final int minUsableIndelScoreToUse, final byte flatDeletionPenalty, final byte flatInsertionPenalty) {
        this.fbargs = flowBasedArgs;
        this.log10globalReadMismappingRate = log10globalReadMismappingRate;
        this.expectedErrorRatePerBase = expectedErrorRatePerBase;
        this.readDisqualificationScale = readDisqualificationScale;
        this.dynamicReadDisqualification = dynamicReadDisqualification;
        this.pcrErrorModel = pcrErrorModel;
        this.flowPairHMM = new FlowBasedPairHMM();
        this.dragstrParams = dragstrParams;
        this.constantGCP = constantGCP;
        this.flatDeletionPenalty = flatDeletionPenalty;
        this.flatInsertionPenalty = flatInsertionPenalty;
        this.minUsableIndelScoreToUse = (byte)minUsableIndelScoreToUse;


        initializePCRErrorModel();
    }

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

        result.normalizeLikelihoods(log10globalReadMismappingRate, true);
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
            final double maxErrorsForRead = Math.max(3.0, Math.ceil(read.getLength() * expectedErrorRate));
            final double maxCatastrophicErrorsForRead = Math.max(2.0, Math.ceil(read.getLength() * catastrophicErrorRate));
            return maxErrorsForRead * log10ErrorRate + maxCatastrophicErrorsForRead*catastrophicErrorRate;
        };
    }

    private byte[] pcrIndelErrorModelCache;

    private void initializePCRErrorModel() {

        inputScoreImputator = dragstrParams == null
                ? NonSymmetricalPairHMMInputScoreImputator.newInstance(constantGCP, flatInsertionPenalty, flatDeletionPenalty)
                : DragstrPairHMMInputScoreImputator.of(dragstrParams) ;

        if ( !pcrErrorModel.hasRateFactor() ) {
            return;
        }

        pcrIndelErrorModelCache = new byte[ReadLikelihoodCalculationEngine.MAX_REPEAT_LENGTH + 1];

        final double rateFactor = pcrErrorModel.getRateFactor();

        for( int i = 0; i <= ReadLikelihoodCalculationEngine.MAX_REPEAT_LENGTH; i++ ) {
            pcrIndelErrorModelCache[i] = getErrorModelAdjustedQual(i, rateFactor, minUsableIndelScoreToUse);
        }

    }


    // TODO these methods are ripped whole cloth from the PairHMMLikelihoodsCalculationEngine and should be uinfied/fixed if possible
    @VisibleForTesting
    void applyPCRErrorModel( final byte[] readBases, final byte[] readInsQuals, final byte[] readDelQuals ) {
        if ( pcrErrorModel == PairHMMLikelihoodCalculationEngine.PCRErrorModel.NONE ) {
            return;
        }

        for ( int i = 1; i < readBases.length; i++ ) {
            final int repeatLength = ReadLikelihoodCalculationEngine.findTandemRepeatUnits(readBases, i-1).getRight();
            readInsQuals[i-1] = (byte) Math.min(0xff & readInsQuals[i - 1], 0xff & pcrIndelErrorModelCache[repeatLength]);
            readDelQuals[i-1] = (byte) Math.min(0xff & readDelQuals[i - 1], 0xff & pcrIndelErrorModelCache[repeatLength]);
        }
    }

    private static void capMinimumReadIndelQualities(final byte[] readInsQuals, final byte[] readDelQuals, final byte minUsableQualScore) {
        for( int i = 0; i < readInsQuals.length; i++ ) {
            readInsQuals[i] = setToFixedValueIfTooLow( readInsQuals[i], minUsableQualScore,       minUsableQualScore );
            readDelQuals[i] = setToFixedValueIfTooLow( readDelQuals[i], minUsableQualScore,       minUsableQualScore );
        }
    }

    private static byte setToFixedValueIfTooLow(final byte currentVal, final byte minQual, final byte fixedQual){
        return currentVal < minQual ? fixedQual : currentVal;
    }


    static byte getErrorModelAdjustedQual(final int repeatLength, final double rateFactor, final byte minUsableIndelScoreToUse) {
        return (byte) Math.max(minUsableIndelScoreToUse, MathUtils.fastRound(INITIAL_QSCORE - Math.exp(repeatLength / (rateFactor * Math.PI)) + 1.0));
    }

    /**
     * Initialize our flowPairHMM with parameters appropriate to the haplotypes and reads we're going to evaluate
     *
     * After calling this routine the PairHMM will be configured to best evaluate all reads in the samples
     * against the set of haplotypes
     *  @param haplotypes a non-null list of haplotypes
     * @param perSampleReadList a mapping from sample -> reads
     */
    private void initializeFlowPairHMM(final List<FlowBasedHaplotype> haplotypes, final List<FlowBasedRead> perSampleReadList) {
        initializePCRErrorModel();

        final int readMaxLength = perSampleReadList.stream().mapToInt(FlowBasedRead::getKeyLength).max().orElse(0);
        final int haplotypeMaxLength = haplotypes.stream().mapToInt(h -> h.getKeyLength()).max().orElse(0);

        // initialize arrays to hold the probabilities of being in the match, insertion and deletion cases
        flowPairHMM.initialize(readMaxLength, haplotypeMaxLength);;
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
                max_class = FlowBasedRead.MAX_CLASS;
            } else {
                max_class = Integer.parseInt(mc_string);
            }

            //get flow order for conversion.
            fo = hdr.getReadGroup(rd.getReadGroup()).getFlowOrder();
            if ( fo == null ) {
                throw new GATKException("Unable to perform flow based alignment without the flow order information");
            }
            if (fbargs.flowOrderCycleLength != 4) {
                throw new GATKException("FlowBasedHMMEngine requires flow order of 4 elements but cycle length is specified as " + fbargs.flowOrderCycleLength);
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

            //TODO This imputation code will eventually need to be turned into something real based on dragstr:
            final PairHMMInputScoreImputation inputScoreImputation = inputScoreImputator.impute(tmp);
            final byte[] readInsQuals = inputScoreImputation.insOpenPenalties();
            final byte[] readDelQuals = inputScoreImputation.delOpenPenalties();
            final byte[] overallGCP = inputScoreImputation.gapContinuationPenalties();

            applyPCRErrorModel(tmp.getBases(), readInsQuals, readDelQuals);
            capMinimumReadIndelQualities(readInsQuals, readDelQuals, minUsableIndelScoreToUse);

            tmp.setReadInsQuals(readInsQuals);
            tmp.setReadDelQuals(readDelQuals);
            tmp.setOverallGCP(overallGCP);

            processedReads.add(tmp);
        }

        //same for the haplotypes - each haplotype is converted to FlowBasedHaplotype
        FlowBasedHaplotype fbh;

        flowOrder = FlowBasedAlignmentEngine.findFlowOrderForReadGroups(hdr, flowOrder, fbargs);

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
            fbr.applyBaseClipping(Math.max(0, diffLeft), Math.max(diffRight, 0), false);
        }
        initializeFlowPairHMM(processedHaplotypes, processedReads);


        flowPairHMM.computeLog10LikelihoodsFlowBased(likelihoods, processedReads, processedHaplotypes);
    }

    @Override
    public void close() {}

}