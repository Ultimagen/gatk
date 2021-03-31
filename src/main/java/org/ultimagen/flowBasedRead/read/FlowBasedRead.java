package org.ultimagen.flowBasedRead.read;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;
import org.apache.commons.lang.ArrayUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.Tail;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.clipping.ClippingOp;
import org.broadinstitute.hellbender.utils.clipping.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.ultimagen.flowBasedRead.utils.Direction;
import org.ultimagen.flowBasedRead.utils.FlowBasedAlignmentArgumentCollection;

import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
/*
Adds flow information to the usual GATKRead. In addition to the usual read data this class keeps flowMatrix,
that contains probabilities for alternative hmer calls.

Main function deals with parsing flow-specific QUAL representation readBaseMatrixProb.
Note that there is a lot of code that deals with other varous formats of the representation (e.g. when the matrix
is coded in the tags of the BAM and is given in flow space). This code is not used in production, but was used in
development and testing
*/

public class FlowBasedRead extends SAMRecordToGATKReadAdapter implements GATKRead, FlowBasedReadInterface, Serializable {

    private static final long serialVersionUID = 42L;
    GATKRead read = null;
    private static int N_ASCII=78;

    private SAMRecord samRecord;
    private byte[] forwardSequence;
    private byte[] key;
    private int [] flow2base;
    private int maxHmer;
    private byte[] flowOrder;
    private double[][] flowMatrix;
    private boolean validKey;
    private Direction direction = Direction.SYNTHESIS;
    private boolean trimmedToHaplotype = false;
    private int trimLeftBase = 0 ;
    private int trimRightBase = 0 ;
    static private final int MINIMAL_READ_LENGTH = 10; // check if this is the right number
    static private final int MAXIMAL_MAXHMER = 100;
    private final double MINIMAL_CALL_PROB = 0.1;
    private final FlowBasedAlignmentArgumentCollection fbargs;
    private final Logger logger = LogManager.getLogger(this.getClass());
    static private int[] flowMatrixModsInstructions = new int[MAXIMAL_MAXHMER];

    /**
     * Creating flow based read when the argument collection is not defined.
     * @param samRecord Record
     * @param _flowOrder flow order stirng (one cycle)
     * @param _maxHmer maximal hmer to store in the flow matrix
     */
    public FlowBasedRead(final SAMRecord samRecord, final String _flowOrder, final int _maxHmer) {
        this(samRecord, _flowOrder, _maxHmer, new FlowBasedAlignmentArgumentCollection());
    }


    /**
     * copy constructor
     * @param other read
     */
    public FlowBasedRead(final FlowBasedRead other) {
        super(other.samRecord);
        forwardSequence = other.forwardSequence.clone();
        key = other.key.clone();
        flow2base = other.flow2base.clone();
        flowOrder = other.flowOrder.clone();
        flowMatrix = other.flowMatrix.clone();
        validKey = other.validKey;
        direction = other.direction;
        maxHmer = other.maxHmer;
        fbargs = other.fbargs;
    }

    /**
     * Constructor from GATKRead. flow order, hmer and arguments
     * @param read GATK read
     * @param flowOrder flow order string (one cycle)
     * @param _maxHmer maximal hmer to keep in the flow matrix
     * @param fbargs arguments that control resolution etc. of the flow matrix
     */
    public FlowBasedRead(final GATKRead read, final String flowOrder, final int _maxHmer, final FlowBasedAlignmentArgumentCollection fbargs) {
        this(read.convertToSAMRecord(null), flowOrder, _maxHmer, fbargs);
        this.read = read;

    }

    /**
     * Same as above but constructs from SAMRecord
     * @param samRecord record from SAM file
     * @param _flowOrder flow order (single cycle)
     * @param _maxHmer maximal hmer to keep in the flow matrix
     * @param fbargs arguments that control resoltion of the flow matrix
     */
    public FlowBasedRead(final SAMRecord samRecord, final String _flowOrder, final int _maxHmer, final FlowBasedAlignmentArgumentCollection fbargs) {
        super(samRecord);
        this.fbargs = fbargs;
        maxHmer = _maxHmer;
        this.samRecord = samRecord;
        forwardSequence = getForwardSequence();

        //supports old format, where the matrix is stored in flow space in the record
        if ( samRecord.hasAttribute("kr") )
            readFlowMatrix(_flowOrder);
        // supports FASTQ-like format
        else {
            if (samRecord.hasAttribute("ti")) {
                readBaseMatrixRecal(_flowOrder);
            } else if (samRecord.hasAttribute("tp")) {
                readBaseMatrixProb(_flowOrder);
            }
        }


        //Spread boundary flow probabilities when the read is unclipped
        //in this case the value of the hmer is uncertain
        if (CigarUtils.countClippedBases(samRecord.getCigar(), Tail.LEFT, CigarOperator.HARD_CLIP) == 0){
            _spreadFlowProbs(findFirstNonZero(key));
        }
        if (CigarUtils.countClippedBases(samRecord.getCigar(), Tail.RIGHT, CigarOperator.HARD_CLIP) == 0){
            _spreadFlowProbs(findLastNonZero(key));
        }


        if ( logger.isDebugEnabled() ) {
            logger.debug("cons: name: " + samRecord.getReadName()
                    + " len: " + samRecord.getReadLength()
                    + " loc: " + samRecord.getStart() + "-" + samRecord.getEnd()
                    + " rev: " + isReverseStrand()
                    + " cigar:" + samRecord.getCigarString());
            logger.debug("     bases: " + new String(samRecord.getReadBases()));
            logger.debug("       key: " + keyAsString(key));
        }

        validateSequence();
    }

    //since the last unclipped flow is uncertain (we give high probabilities to
    //also hmers higher than the called hmer)
    private void _spreadFlowProbs(final int flowToSpread) {
        if (flowToSpread<0) //boundary case when all the key is zero
            return;

        final int call = key[flowToSpread];
        if (call==0){
            throw new GATKException.ShouldNeverReachHereException("Boundary key value should not be zero for the spreading");
        }

        final int numberToFill = maxHmer - call+1;
        double total = 0;
        for (int i = call; i < maxHmer+1; i++)
            total += flowMatrix[i][flowToSpread];
        final double fillProb = Math.max(total / numberToFill, fbargs.fillingValue);
        for (int i = call; i < maxHmer+1; i++){
            flowMatrix[i][flowToSpread] = fillProb;
        }
    }



    // convert qualities and ti tag to flow matrix
    private void readBaseMatrixRecal(final String _flowOrder) {

       // generate key (base to flow space)
        setDirection(Direction.REFERENCE);  // base is always in reference/alignment direction
        key = base2key(samRecord.getReadBases(), _flowOrder, 1000);
        getKey2Base();
        flowOrder = getFlow2Base(_flowOrder, key.length);

       // initialize matrix
        flowMatrix = new double[maxHmer+1][key.length];
        for (int i = 0 ; i < maxHmer+1; i++) {
            for (int j = 0 ; j < key.length; j++ ){
                flowMatrix[i][j] = fbargs.fillingValue;;
            }
        }

        // access qual, convert to flow representation
        final byte[]      quals = samRecord.getBaseQualities();
        final byte[]      ti = samRecord.getByteArrayAttribute("ti");
        final double[]    probs = new double[quals.length];
        for ( int i = 0 ; i < quals.length ; i++ ) {
            final double q = quals[i];
            final double p = QualityUtils.qualToErrorProb(q);
            probs[i] = p*2;
        }

        // apply key and qual/ti to matrix
        int     qualOfs = 0;
        for ( int i = 0 ; i < key.length ; i++ ) {
            final byte        run = key[i];

            // the probability is not divided by two for hmers of length 1
            if ( run == 1 ) {
                probs[qualOfs] = probs[qualOfs]/2;
            }

            //filling the probability for the called hmer (not reported by the quals
            if ( run <= maxHmer ) {
                flowMatrix[run][i] = (run > 0) ? (1 - probs[qualOfs]) : 1;
                //require a prob. at least 0.1
                flowMatrix[run][i] = Math.max(MINIMAL_CALL_PROB, flowMatrix[run][i]);

            }

            if ( run != 0 ) {
                if ( quals[qualOfs] != 40 ) {
                    final int     run1 = (ti[qualOfs] == 0) ? (run - 1) : (run + 1);
                    if (( run1 <= maxHmer ) && (run <= maxHmer)){
                        flowMatrix[run1][i] = probs[qualOfs] / flowMatrix[run][i];
                    }
                    if (run <= maxHmer) {
                        flowMatrix[run][i] /= flowMatrix[run][i]; // for comparison to the flow space - probabilities are normalized by the key's probability
                    }
                }
                qualOfs += run;
            }

        }

        //this is just for tests of all kinds of
        applyFilteringFlowMatrix();
    }


    //This is the code for parsing the current BAM format (with TP tag)
    private void readBaseMatrixProb(final String _flowOrder) {

        // generate key (base to flow space)
        setDirection(Direction.REFERENCE);  // base is always in reference/alignment direction

        key = base2key(samRecord.getReadBases(), _flowOrder, 1000);
        getKey2Base();
        flowOrder = getFlow2Base(_flowOrder, key.length);

        // initialize matrix
        flowMatrix = new double[maxHmer+1][key.length];
        for (int i = 0 ; i < maxHmer+1; i++) {
            for (int j = 0 ; j < key.length; j++ ){
                flowMatrix[i][j] = fbargs.fillingValue;
            }
        }

        // access qual, convert to ultima representation
        final byte[]      quals = samRecord.getBaseQualities();
        final byte[]      tp = samRecord.getSignedByteArrayAttribute("tp");
        final double[]    probs = new double[quals.length];
        for ( int i = 0 ; i < quals.length ; i++ ) {
            final double q = quals[i];
            final double p = Math.pow(10, -q/10);
            probs[i] = p;
        }

        // apply key and qual/tp to matrix
        int     qualOfs = 0; //converts between base -> flow
        for ( int i = 0 ; i < key.length ; i++ ) {
            final byte        run = key[i];
            if (run > 0) {
                parseSingleHmer(probs, tp, i, run, qualOfs);
            }
            double totalErrorProb = 0;

            for (int k=0; k < maxHmer; k++ ){
                totalErrorProb += flowMatrix[k][i];
            }
            final double callProb = Math.max(MINIMAL_CALL_PROB, 1-totalErrorProb);
            // the probability in the recalibration is not divided by two for hmers of length 1
            flowMatrix[Math.min(run, maxHmer)][i] = callProb;
            qualOfs+=run;
        }
        applyFilteringFlowMatrix();
    }


    //convert qualities from the single hmer to a column in a flow matrix
    private void parseSingleHmer(final double[] probs, final byte[] tp, final int flowIdx,
                                 final byte flowCall, final int qualOfs){
        for (int i = qualOfs ; i < qualOfs+flowCall; i++) {
            if (tp[i]!=0) {
                final int loc = Math.max(Math.min(flowCall+tp[i], maxHmer),0);
                if (flowMatrix[loc][flowIdx] == fbargs.fillingValue) {
                    flowMatrix[loc][flowIdx] = probs[i];
                } else {
                    flowMatrix[loc][flowIdx] += probs[i];
                }
            }
        }
    }

    public String getFlowOrder() {
        return new String(Arrays.copyOfRange(flowOrder, 0, Math.min(fbargs.flowOrderCycleLength,flowOrder.length)));
    }


    public int getMaxHmer() {
        return maxHmer;
    }

    public int getNFlows() {
        return key.length;
    }

    public Direction getDirection(){
        return direction;
    }


    private byte[] getForwardSequence(){
        if (!isReverseStrand()) {
            return samRecord.getReadBases();
        } else {
            final byte[] result = new byte[samRecord.getReadBases().length];
            System.arraycopy(samRecord.getReadBases(), 0, result, 0, result.length);
            SequenceUtil.reverseComplement(result);
            return result;
        }
    }


    private int[] getAttributeAsIntArray(final String attributeName, final boolean isSigned) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        final Object attributeValue = this.samRecord.getAttribute(attributeName);

        if (attributeValue == null) {
            return null;
        } else if (attributeValue instanceof byte[]) {
            final byte[] tmp = (byte[]) attributeValue;
            final int[] ret = new int[tmp.length];
            for (int i = 0; i < ret.length; i++)
                if (!isSigned) ret[i] = tmp[i]&0xff;
                else ret[i]=tmp[i]; //converting signed byte to unsigned
            return Arrays.copyOf(ret, ret.length);
        } else if ((attributeValue instanceof int[])) {
            final int[] ret = (int[]) attributeValue;
            return Arrays.copyOf(ret, ret.length);
        } else if  (attributeValue instanceof short[]) {
            final short [] tmp = (short[]) attributeValue;
            final int[] ret = new int[tmp.length];
            for (int i = 0 ; i < tmp.length; i++ )
                ret[i] = tmp[i];
            return Arrays.copyOf(ret, ret.length);
        }else {
            throw new GATKException.ReadAttributeTypeMismatch(attributeName, "integer array");
        }
    }


    private void validateSequence(){
        for (final byte b : key) {
            if (b > maxHmer - 1) {
                validKey = false;
            }
        }
        validKey = true;
    }

    public boolean isValid() {
        return validKey;
    }

    public double getProb(final int flow, final int hmer) {
        return flowMatrix[Math.min(hmer, maxHmer)][flow];
    }

    // this applies clipping when the flow matrix is in flow space. Does nothing in base space
    public void applyAlignment(){
        if ((getDirection() == Direction.SYNTHESIS) && ( isReverseStrand() )) {
            flipMatrix();
            reverse(key, key.length);
            getKey2Base();
            SequenceUtil.reverseComplement(flowOrder);

        }

        final boolean isBase = isBaseFormat();
        final int[] basePair = {0, 0};
        final int[] clipLeftPair = !isBase ? findLeftClipping() : basePair;
        final int[] clipRightPair = !isBase ? findRightClipping() : basePair;
        final int clipLeft = clipLeftPair[0];
        final int leftHmerClip = clipLeftPair[1];
        final int clipRight = clipRightPair[0];
        final int rightHmerClip = clipRightPair[1];

        applyClipping(clipLeft, leftHmerClip, clipRight, rightHmerClip);

        setDirection(Direction.REFERENCE);

    }

    private boolean isBaseFormat() {
       return samRecord.hasAttribute("ti") || samRecord.hasAttribute("tp");
    }


    private void reverse(final int []a, final int n)
    {
        int i, k, t;
        for (i = 0; i < n / 2; i++) {
            t = a[i];
            a[i] = a[n - i - 1];
            a[n - i - 1] = t;
        }

    }


    private void reverse(final byte []a, final int n)
    {
        int i, k;
        byte t;
        for (i = 0; i < n / 2; i++) {
            t = a[i];
            a[i] = a[n - i - 1];
            a[n - i - 1] = t;
        }

    }

    private void reverse(final double []a, final int n)
    {
        int i, k;
        double t;
        for (i = 0; i < n / 2; i++) {
            t = a[i];
            a[i] = a[n - i - 1];
            a[n - i - 1] = t;
        }

    }


    // code for reading BAM format where the flow matrix is stored in sparse representation in kr,kf,kh and kd tags
    // used for development of the new basecalling, but not in production code
    private void readFlowMatrix(final String _flowOrder) {

        key = getAttributeAsByteArray("kr");

        // creates a translation from flow # to base #
        getKey2Base();

        // create a translation from
        flowOrder = getFlow2Base(_flowOrder, key.length);

        flowMatrix = new double[maxHmer+1][key.length];
        for (int i = 0 ; i < maxHmer+1; i++) {
            for (int j = 0 ; j < key.length; j++ ){
                flowMatrix[i][j] = fbargs.fillingValue;
            }
        }

        byte [] kh = getAttributeAsByteArray( "kh" );
        int [] kf = getAttributeAsIntArray("kf", false);
        int [] kd = getAttributeAsIntArray( "kd", true);

        final byte [] key_kh = key;
        final int [] key_kf = new int[key.length];
        for ( int i = 0 ; i < key_kf.length ; i++)
            key_kf[i] = i;
        final int [] key_kd = new int[key.length];

        kh = ArrayUtils.addAll(kh, key_kh);
        kf = ArrayUtils.addAll(kf, key_kf);
        kd = ArrayUtils.addAll(kd, key_kd);

        quantizeProbs(kd);

        final double [] kdProbs = phredToProb(kd);
        fillFlowMatrix( kh, kf, kdProbs);
        applyFilteringFlowMatrix();
        validateSequence();
    }

    private void fillFlowMatrix(final byte [] kh, final int [] kf,
                                final double [] kdProbs ) {
        for ( int i = 0 ; i < kh.length; i++ ) {
            // normal matrix filling
            final int     pos = kf[i];
            final int     hmer = kh[i] & 0xff;
            if (hmer > maxHmer){
                flowMatrix[maxHmer][pos] = Math.max(flowMatrix[maxHmer][pos], kdProbs[i]);
            } else {
                flowMatrix[hmer][pos] = Math.max(flowMatrix[hmer][pos], kdProbs[i]);
            }

            // mod instruction implementation
            final int hmer2 = flowMatrixModsInstructions[hmer];
            if ( hmer2 != 0 ) {
                flowMatrix[hmer2][pos] = Math.max(flowMatrix[hmer2][pos], kdProbs[i]);

                // if we are copying bacwards, zero out source
                if ( hmer > hmer2 )
                    flowMatrix[hmer][pos] = 0;
            }
        }

    }


    private void flipMatrix() {
        for ( int i = 0 ; i < flowMatrix.length; i++) reverse(flowMatrix[i], flowMatrix[i].length);
    }

    private int findFirstNonZero(final byte[] array){
        int result = -1;
        for (int i = 0 ; i < array.length; i++){
            if (array[i]!=0) {
                result = i;
                break;
            }
        }
        return result;
    }

    private int findLastNonZero(final byte[] array){
        int result = -1;
        for (int i = array.length-1 ; i >= 0; i--){
            if (array[i]!=0) {
                result = i;
                break;
            }
        }
        return result;
    }

    private void shiftColumnUp(final double[][] matrix, final int colnum, final int shift) {
        for (int i = 0; i < matrix.length - shift; i ++ ) {
            matrix[i][colnum] = matrix[i+shift][colnum];
        }
        for (int i = matrix.length - shift; i < matrix.length; i ++ ) {
            matrix[i][colnum] = 0;
        }

    }

    public void setDirection(final Direction dir ) {
        direction = dir;
    }

    /**
     * Converts base space sequence to flow space
     * @param bases base space sequence
     * @param flowOrder flow order
     * @param clipping maximal flow value to output (at most 127).
     * @return Array of flow values
     */
    static protected byte[] base2key(final byte[] bases, final String flowOrder, int clipping) {

        clipping = Math.min(clipping, 127);

        final int[]       intKeys = base2key(bases, flowOrder);
        final byte[]      byteKeys = new byte[intKeys.length];
        int         i = 0;

        for ( final int intKey : intKeys ) {
            byteKeys[i++] = (byte)((intKey < clipping) ? intKey : clipping);
        }

        return byteKeys;
    }

    /**
     * Converts base space sequence to flow space
     * @param bases base space sequence
     * @param flowOrder flow order
     * @return Array of flow values
     */

    static protected int[] base2key(final byte[] bases, final String flowOrder){

        final ArrayList<Integer> result = new ArrayList<>();
        final byte[] flowOrderBytes = flowOrder.getBytes();
        int loc = 0;
        int flowNumber = 0 ;
        final int period = flowOrderBytes.length;
        while ( loc < bases.length ) {
            final byte flowBase = flowOrderBytes[flowNumber%period];
            if ((bases[loc]!=flowBase) && ( bases[loc]!=N_ASCII)) {
                result.add(0);
            } else {
                int count = 0;
                while ( ( loc < bases.length) && ((bases[loc]==flowBase) || (bases[loc]==N_ASCII)) ){
                    loc++;
                    count ++;
                }
                result.add(count);

            }
            flowNumber++;
        }
        final int[] ret = new int[result.size()];
        for (int i = 0; i < result.size(); i++) {
            ret[i] = result.get(i);
        }
        return ret;
    }

    // For every flow of the key output the index of the last base that was output prior to this flow
    private void getKey2Base() {
        flow2base = new int[key.length];
        flow2base[0] = -1;
        for (int i = 1 ; i < flow2base.length; i++) {
            flow2base[i] = flow2base[i-1] + key[i-1];
        }
    }

    /**
     * For every flow of the key output the index of the last base that was output prior to this flow
     * @param key
     * @return array
     */
    static protected int[] getKey2Base(final byte[] key){
        final int[] result = new int[key.length];
        result[0] = -1;
        for (int i = 1 ; i < result.length; i++) {
            result[i] = result[i-1] + key[i-1];
        }
        return result;
    }

    /**
     * For every flow of the key output the index of the last base that was output prior to this flow
     * @param key
     * @return array
     */
    static protected int[] getKey2Base(final int[] key) {
        final int[] result = new int[key.length];
        result[0] = -1;
        for (int i = 1; i < result.length; i++) {
            result[i] = result[i - 1] + key[i - 1];
        }
        return result;

    }

    /**
     * For every flow of the key output the nucleotide that is being read for this flow
     * @param flowOrder
     * @param expectedLength the length of the key (key is not provided)
     * @return array of bases
     */

    static protected byte[] getFlow2Base(final String flowOrder, final int expectedLength) {
        final byte[] result = new byte[expectedLength] ;
        for ( int i = 0; i < result.length; i++ ) {
            result[i] = (byte)flowOrder.charAt(i%flowOrder.length());
        }
        return result;
    }

    //trims base-spaced reads. Usually not needed, but kept for completeness
    public void applyBaseClipping(final int clipLeftBase, final int clipRightBase){
        final int[] clipLeftPair = findLeftClipping(clipLeftBase);
        final int[] clipRightPair = findRightClipping(clipRightBase);
        final int clipLeft = clipLeftPair[0];
        final int leftHmerClip = clipLeftPair[1];
        final int clipRight = clipRightPair[0];
        final int rightHmerClip = clipRightPair[1];
        if (getLength() - clipLeftBase - clipRightBase < MINIMAL_READ_LENGTH) {
            trimmedToHaplotype = true;
            validKey =false;
            trimLeftBase =clipLeftBase;
            trimRightBase = clipRightBase;
        } else {
            applyClipping(clipLeft, leftHmerClip, clipRight, rightHmerClip);
            trimmedToHaplotype = true;
            trimLeftBase = clipLeftBase;
            trimRightBase = clipRightBase;
        }
    }

    private void applyClipping(int clipLeft, final int leftHmerClip, int clipRight, final int rightHmerClip){
        if ((clipLeft < 0) || (clipRight < 0)  || (clipLeft >= getKeyLength() ) || ( clipRight >= getKeyLength())) {
            throw new GATKException.ShouldNeverReachHereException("Weird read clip calculated");
            //return 1;
        }

        if ((leftHmerClip < 0) || (rightHmerClip < 0)  || (leftHmerClip >= 14 ) || ( rightHmerClip >= 14)) {
            throw new GATKException.ShouldNeverReachHereException("Weird read clip calculated");
            //return 1;
        }

        final int originalLength = key.length;

        key[clipLeft]-=leftHmerClip;
        boolean shiftLeft = true;
        if ( (clipLeft >= 0) || ( leftHmerClip >= 0 )  ) {
            while (key[clipLeft] == 0) {
                clipLeft += 1 ;
                shiftLeft = false;
            }
        }
        key[key.length - clipRight-1] -= rightHmerClip;
        boolean shiftRight = true;
        if ( (clipRight >= 0) || ( rightHmerClip >= 0 )  ) {
            while (key[originalLength - 1- clipRight] == 0) {
                clipRight += 1 ;
                shiftRight = false;
            }
        }

        key = Arrays.copyOfRange(key, clipLeft, originalLength - clipRight);
        getKey2Base();
        flowOrder = Arrays.copyOfRange(flowOrder, clipLeft, originalLength - clipRight);

        final double [][] newFlowMatrix = new double[flowMatrix.length][originalLength - clipLeft - clipRight] ;
        for ( int i = 0 ; i < newFlowMatrix.length; i++) {
            newFlowMatrix[i] = Arrays.copyOfRange(flowMatrix[i], clipLeft, originalLength - clipRight);
        }

        flowMatrix = newFlowMatrix;
        if (shiftLeft) {
            shiftColumnUp(flowMatrix, 0, leftHmerClip);
        }

        if (shiftRight) {
            shiftColumnUp(flowMatrix, flowMatrix[0].length-1, rightHmerClip);
        }

        //Spread boundary flow probabilities for the boundary hmers of the read
        //in this case the value of the genome hmer is uncertain
        _spreadFlowProbs(findFirstNonZero(key));
        _spreadFlowProbs(findLastNonZero(key));
    }

    private int[] findLeftClipping() {
        final List<CigarElement> cigar = getCigarElements();
        final int[] result = new int[2];
        if (cigar.size() == 0 ) {
            return result;
        }

        final CigarElement start = cigar.get(0);
        if (start.getOperator() != CigarOperator.H) {
            return result;
        }

        final int basesClipped = start.getLength();
        return findLeftClipping(basesClipped);
    }


    private int[] findLeftClipping(final int basesClipped){
        final int[] result = new int[2];
        if (basesClipped==0){
            return result;
        }

        int stopClip = 0;
        for (int i = 0 ; i < flow2base.length; i++ ) {

            if (flow2base[i] + key[i] >= basesClipped) {
                stopClip = i;
                break;
            }
        }
        final int hmerClipped = basesClipped - flow2base[stopClip] - 1;
        result[0] = stopClip;
        result[1] = hmerClipped;
        return result;
    }

    private int[] findRightClipping() {
        final List<CigarElement> cigar = getCigarElements();
        final int[] result = new int[2];
        if (cigar.size() == 0 ) {
            result[0] = 0;
            result[1] = 0;
            return result;
        }

        final CigarElement end = cigar.get(cigar.size()-1);
        if (end.getOperator() != CigarOperator.H) {
            result[0] = 0;
            result[1] = 0;
            return result;
        }

        final int basesClipped = end.getLength();

        return findRightClipping(basesClipped);
    }

    private int[] findRightClipping(final int basesClipped) {
        final int[] result = new int[2];
        if (basesClipped==0){
            return result;
        }

        final int index =0 ;
        int stopClip = 0;

        final byte[] rkey = new byte[key.length];
        for (int i = 0 ; i < key.length; i++ ){
            rkey[i] = key[key.length-1-i];
        }

        final int[] rflow2base = getKey2Base(rkey);
        for (int i = 0 ; i < rflow2base.length; i++ ) {
            if (rflow2base[i] + rkey[i] >= basesClipped) {
                stopClip = i;
                break;
            }
        }

        final int hmerClipped = basesClipped - rflow2base[stopClip] - 1;
        result[0] = stopClip;
        result[1] = hmerClipped;
        return result;
    }


    public void writeKey(final FileWriter oos)
            throws IOException {
        for (int i = 0; i < key.length; i++)
            oos.write(key[i]+"\n");

    }

    /**
     * Flow matrix logger
     * @param oos
     * @throws IOException
     */
    public void writeMatrix(final OutputStreamWriter oos)
            throws IOException {
        final DecimalFormat formatter = new DecimalFormat("0.0000", DecimalFormatSymbols.getInstance(Locale.ENGLISH));

        final byte[]      bases = samRecord.getReadBases();
        int         basesOfs = 0;
        final byte[]      quals = samRecord.getBaseQualities();
        final byte[]      ti = samRecord.hasAttribute("ti") ? samRecord.getByteArrayAttribute("ti") : (new byte[key.length*3]);
        if ( isReverseStrand() )
        {
            reverse(quals, quals.length);
            reverse(ti, ti.length);
        }

        for ( int col = 0 ; col < key.length ; col++ ) {
            oos.write("C,R,F,B,Bi,Q,ti\n");
            final byte base = (key[col] != 0) ? (basesOfs < bases.length ? bases[basesOfs] : (byte)'?') : (byte)'.';
            final String bi = (key[col] != 0) ? Integer.toString(basesOfs) : ".";
            final String q = (key[col] != 0) ? Integer.toString(quals[basesOfs]) : ".";
            final String Ti = (key[col] != 0) ? Integer.toString(ti[basesOfs]) : ".";
            for (int row = 0; row < flowMatrix.length; row++) {
                final String s = formatter.format(flowMatrix[row][col]);
                oos.write(String.format("%d,%d,%d,%c,%s,%s,%s,%s %s\n", col, row, key[col], base, bi, q, Ti, isReverseStrand() ? "r" : ".", s));
            }
            if ( key[col] != 0 )
                basesOfs +=  key[col];
            oos.write("\n");
        }

    }

    /**
     * Prints the key as character-encoded string
     * @param bytes  key array
     * @return encoding
     */
    static public String keyAsString(final byte[] bytes)
    {
        final StringBuilder   sb = new StringBuilder();

        for ( final byte b : bytes )
            sb.append((char)((b < 10) ? ('0' + b) : ('A' + b - 10)));

        return sb.toString();
    }

    /**
     * Prints the key as character-encoded string
     * @param ints (key array)
     * @return encoded string
     */
    static public String keyAsString(final int[] ints)
    {
        final StringBuilder   sb = new StringBuilder();

        for ( final int i : ints )
            sb.append((char)((i < 10) ? ('0' + i) : ('A' + i - 10)));

        return sb.toString();
    }


    public byte [] getFlowOrderArray() {
        return flowOrder;
    }

    public int getKeyLength() {
        return key.length;
    }

    public byte[] getKey() {
        return key;
    }

    /**
     * Number of total bases that the flow based key generates
     * @return number of bases
     */
    public int totalKeyBases()  {
        int sum = 0 ;
        for (int i = 0 ; i < key.length; i++){
            sum += key[i];
        }
        return sum;
    }

    public int seqLength(){
        return forwardSequence.length;
    }
    public boolean isTrimmedToHaplotype() {
        return trimmedToHaplotype;
    }

    public int getTrimmedStart() {
        return trimLeftBase + getStart();
    }
    public int getTrimmedEnd() {
        return getEnd() - trimRightBase;
    }

    //functions that take care of simulating base format
    //they perform modifications on the flow matrix that are defined in applyFilteringFlowMatrix

    //this function was only applied when we tested what is the necessary information to be reported in the flow matrix
    private void applyFilteringFlowMatrix(){

        if (fbargs.disallow_larger_probs) {
            removeLargeProbs();
        }

        if (fbargs.remove_longer_than_one_indels) {
            removeLongIndels( key );
        }

        if (fbargs.remove_one_to_zero_probs) {
            removeOneToZeroProbs(key);
        }

        if ((fbargs.lump_probs)) {
            lumpProbs();
        }
        clipProbs();

        if (fbargs.symmetric_indels) {
            smoothIndels(key);
        }
        if (fbargs.only_ins_or_del) {
            reportInsOrDel(key);
        }

        if ((fbargs.retainMaxNProbs)){
            reportMaxNProbsHmer(key);
        }

    }

    private void clipProbs() {
        for ( int i = 0 ; i < getMaxHmer(); i++ ) {
            for ( int j =0; j < getNFlows(); j++) {
                if ((flowMatrix[i][j] < fbargs.probability_ratio_threshold) &&
                        (key[j]!=i)) {
                    flowMatrix[i][j] = fbargs.fillingValue;
                }
            }
        }
    }

    private void removeLargeProbs(){
        for (int i = 0; i < getNFlows(); i++){
            for (int j = 0 ; j < getMaxHmer()+1; j++) {
                if (flowMatrix[j][i] > 1) {
                    flowMatrix[j][i] = 1;
                }
            }
        }
    }

    private void removeLongIndels(final byte [] key_kh ){
        for ( int i = 0 ; i < getNFlows(); i++ ) {
            for (int j = 0; j < getMaxHmer()+1; j++){
                if (Math.abs(j-key_kh[i])>1){
                    flowMatrix[j][i] = fbargs.fillingValue;
                }
            }
        }
    }

    private void removeOneToZeroProbs(final byte [] key_kh) {
        for (int i = 0 ; i < getNFlows(); i++){
            if (key_kh[i] == 0){
                for (int j = 1; j < getMaxHmer()+1; j++){
                    flowMatrix[j][i]=fbargs.fillingValue;
                }
            }
        }
    }



    private void quantizeProbs(final int [] kd_probs ) {
        final int nQuants = fbargs.probability_quantization;
        final double bin_size = 6*fbargs.probabilityScalingFactor/(float)nQuants;
        for ( int i = 0 ; i < kd_probs.length; i++) {
            if (kd_probs[i] <=0)
                continue;
            else {
                kd_probs[i] = (byte)(bin_size * (byte)(kd_probs[i]/bin_size)+1);
            }
        }
    }

    private void quantizeProbs() {

        final int nQuants = fbargs.probability_quantization;
        final double bin_size = 6*fbargs.probabilityScalingFactor/(float)nQuants;
        for (int i = 0 ; i < getNFlows(); i++) {
            for (int j = 0 ; j < getMaxHmer()+1; j ++){
                if ( flowMatrix[j][i] == fbargs.fillingValue)
                    continue;
                if (flowMatrix[j][i]>=1){
                    continue;
                }

                final double origQual = -fbargs.probabilityScalingFactor*Math.log10(flowMatrix[j][i]);
                final byte binnedQual = (byte)(bin_size * (byte)(origQual/bin_size)+1);
                flowMatrix[j][i] = Math.max(Math.pow(10, ((double)binnedQual)/fbargs.probabilityScalingFactor), fbargs.fillingValue);
            }
        }
    }


    private void smoothIndels(final byte [] kr ) {
        for ( int i = 0 ; i < kr.length; i++ ){
            final byte idx = kr[i];
            if (( idx > 1 ) && ( idx < maxHmer) ) {
                final double tmp = (flowMatrix[idx - 1][i] + flowMatrix[idx + 1][i]) / 2;
                flowMatrix[idx - 1][i] = tmp;
                flowMatrix[idx + 1][i] = tmp;
            }
        }
    }

    private void reportInsOrDel(final byte [] kr ) {
        for ( int i = 0 ; i < kr.length; i++ ){
            final byte idx = kr[i];
            if (( idx > 1 ) && ( idx < maxHmer) ) {
                if ((flowMatrix[idx-1][i] > fbargs.fillingValue) && (flowMatrix[idx+1][i] > fbargs.fillingValue)) {
                    final int fixCell = flowMatrix[idx-1][i] > flowMatrix[idx+1][i] ? idx+1 : idx-1;
                    final int otherCell = flowMatrix[idx-1][i] > flowMatrix[idx+1][i] ? idx-1 : idx+1;
                    flowMatrix[fixCell][i] = fbargs.fillingValue;
                }
            }
        }
    }

    private void lumpProbs() {

        for (int i = 0; i < getMaxHmer(); i++) {
            for (int j = 0 ; j < getNFlows(); j ++ ) {
                final int fkey = key[j];
                if (flowMatrix[i][j]<=fbargs.fillingValue) {
                    continue;
                } else {
                    if ( (i - fkey) < -1 ){
                        flowMatrix[fkey-1][j]+=flowMatrix[i][j];
                        flowMatrix[i][j] = fbargs.fillingValue;
                    } else if ((i-fkey) > 1) {
                        flowMatrix[fkey+1][j]+=flowMatrix[i][j];
                        flowMatrix[i][j] = fbargs.fillingValue;
                    }

                }

            }
        }

    }


    private void reportMaxNProbsHmer(final byte [] key) {
        final double [] tmpContainer = new double[maxHmer];
        for (int i = 0 ; i < key.length;i++){

            for (int j = 0 ; j < tmpContainer.length; j++) {
                tmpContainer[j] = flowMatrix[j][i];
            }
            final int k = (key[i]+1)/2;
            final double kth_highest = findKthLargest(tmpContainer, k+1);
            for (int j = 0 ; j < maxHmer; j++)
                if (flowMatrix[j][i] < kth_highest)
                    flowMatrix[j][i] = fbargs.fillingValue;
        }

    }


    private double findKthLargest(final double[] nums, final int k) {
        final PriorityQueue<Double> q = new PriorityQueue<Double>(k);
        for(final double i: nums){
            q.offer(i);

            if(q.size()>k){
                q.poll();
            }
        }

        return q.peek();
    }


    public static void setFlowMatrixMods(final String list) {
        if ( list != null )
        {
            final String[]    toks = list.split(",");
            for ( int i = 0 ; i < toks.length - 1 ; i += 2 )
                flowMatrixModsInstructions[Integer.parseInt(toks[i])] = Integer.parseInt(toks[i+1]);
        }
    }

    private double[] phredToProb(final int [] kq) {
        final double [] result = new double[kq.length];
        for (int i = 0 ; i < kq.length; i++ ) {
            //disallow probabilities below fillingValue
            result[i] = Math.max(Math.pow(10, ((double)-kq[i])/fbargs.probabilityScalingFactor), fbargs.fillingValue);
        }
        return result;
    }

    /**
     * Hard clips uncertain flows (currently four first flows read that often generate noisy base calls
     * @param inputRead GATKREad
     * @param flowOrder flow order string from the SAM header
     * @param fbargs arguments
     * @return read with flowNumUncertainFlows trimmed
     */
    public static GATKRead hardClipUncertainBases(final GATKRead inputRead, final String flowOrder,
                                                  final FlowBasedAlignmentArgumentCollection fbargs ){
        Utils.validateArg(fbargs.flowFirstUncertainFlowBase.length()==1, "First uncertain flow base should be of length 1");
        ReadClipper clipper = new ReadClipper(inputRead);
        final String adjustedFlowOrder = adjustFlowOrderToUncertainFlow(flowOrder, fbargs.flowFirstUncertainFlowBase.charAt(0), fbargs.flowOrderCycleLength);
        if (inputRead.isReverseStrand()) {
            final int nUncertain = nUncertainBases(inputRead, adjustedFlowOrder, fbargs.flowNumUncertainFlows, false);
            clipper.addOp(new ClippingOp(inputRead.getLength()-nUncertain, inputRead.getLength()-1));
        }
        else  {
            final int nUncertain = nUncertainBases(inputRead, adjustedFlowOrder, fbargs.flowNumUncertainFlows, true);
            clipper.addOp(new ClippingOp(0, nUncertain-1));
        }

        return clipper.clipRead(ClippingRepresentation.HARDCLIP_BASES);

    }

    /**
     * Hard clips uncertain flows (currently four first flows read that often generate noisy base calls
     * @param inputRead GATKREad
     * @param samHeader  sam file header (to extract to flow order from
     * @param fbargs arguments
     * @return read with flowNumUncertainFlows trimmed
     */
    public static GATKRead hardClipUncertainBases(final GATKRead inputRead, final SAMFileHeader samHeader,
                                                  final FlowBasedAlignmentArgumentCollection fbargs ){
        ReadClipper clipper = new ReadClipper(inputRead);
        String flowOrder = samHeader.getReadGroup(inputRead.getReadGroup()).getFlowOrder();
        if (flowOrder==null){
            throw new GATKException("Unable to trim uncertain bases without flow order information");
        }
        flowOrder = flowOrder.substring(0,fbargs.flowOrderCycleLength);
        return hardClipUncertainBases(inputRead, flowOrder, fbargs);
    }

    /**
     * Checks if the read has FO tag and thus is flow based
     * @param samHeader header
     * @param inputRead read
     * @return boolean
     */
    public static boolean isFlowBasedData(final SAMFileHeader samHeader, final GATKRead inputRead) {
        String flowOrder = samHeader.getReadGroup(inputRead.getReadGroup()).getFlowOrder();
        if (flowOrder == null) {
            return false;
        }
        return true;
    }

    private static String adjustFlowOrderToUncertainFlow(final String flowOrder,
                                                         final char firstUncertainFlowBase,
                                                         final int flowOrderLength){
        String result = flowOrder + flowOrder;
        final int adjustedStartPos = result.indexOf(firstUncertainFlowBase);
        return result.substring(adjustedStartPos, adjustedStartPos+flowOrderLength);
    }

    //find how many bases are output from uncertain flows
    private static int nUncertainBases(final GATKRead inputRead, final String flowOrder,
                                       final int nUncertainFlows, final boolean isForward){
        byte [] bases;
        if (isForward){
            bases = inputRead.getBases();
        } else {
            bases = ReadUtils.getBasesReverseComplement(inputRead).getBytes();
        }

        final int[] key = base2key(bases, flowOrder);
        final int nTrimFlows = Math.min(nUncertainFlows, key.length);
        int result = 0;
        for (int i = 0 ; i < nTrimFlows; i++){
            result += key[i];
        }
        return result;
    }
}


