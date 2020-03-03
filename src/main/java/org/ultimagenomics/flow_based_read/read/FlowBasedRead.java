package org.ultimagenomics.flow_based_read.read;

import htsjdk.samtools.*;
import org.apache.commons.lang.ArrayUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ultimagenomics.flow_based_read.utils.Direction;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.ultimagenomics.flow_based_read.utils.FlowBasedAlignmentArgumentCollection;

import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Locale;

public class FlowBasedRead extends SAMRecordToGATKReadAdapter implements GATKRead, FlowBasedReadInterface, Serializable {

    private static final long serialVersionUID = 42L;
    private GATKRead read = null;
    private SAMRecord samRecord;
    private byte[] forward_sequence;
    private byte[] key;
    private int [] flow2base;
    private int max_hmer;
    private byte[] flow_order;
    private double[][] flow_matrix;
    private boolean valid_key;
    private Direction direction = Direction.SYNTHESIS;
    private boolean trimmed_to_haplotype = false;
    private int trim_left_base = 0 ;
    private int trim_right_base = 0 ;
    private final int MINIMAL_READ_LENGTH = 10; // check if this is the right number

    private final double ERROR_PROB=1e-4;
    private final double MINIMAL_CALL_PROB = 0.1;
    private final FlowBasedAlignmentArgumentCollection fbargs;
    private final Logger logger = LogManager.getLogger(this.getClass());

    public FlowBasedRead(SAMRecord samRecord, String _flow_order, int _max_hmer) {
        this(samRecord, _flow_order, _max_hmer, new FlowBasedAlignmentArgumentCollection());
    }

    public FlowBasedRead(SAMRecord samRecord, String _flow_order, int _max_hmer, FlowBasedAlignmentArgumentCollection fbargs) {
        super(samRecord);
        this.fbargs = fbargs;
        max_hmer = _max_hmer;
        this.samRecord = samRecord;
        forward_sequence = getForwardSequence();

        if ( samRecord.hasAttribute("kr") )
            readFlowMatrix(_flow_order);
        else
            readBaseMatrix(_flow_order);

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

    static public String keyAsString(byte[] bytes)
    {
        StringBuilder   sb = new StringBuilder();

        for ( byte b : bytes )
            sb.append((char)((b < 10) ? ('0' + b) : ('A' + b - 10)));

        return sb.toString();
    }

    private void readBaseMatrix(String _flow_order) {

       // generate key (base to flow space)
        key = FlowBasedHaplotype.base2key(samRecord.getReadBases(), _flow_order, 1000);
        if ( isReverseStrand() )
            reverse(key, key.length);
        getKey2Base();
        flow_order = getFlow2Base(_flow_order, key.length);
        if ( isReverseStrand() )
            SequenceUtil.reverseComplement(flow_order);

       // initialize matrix
        flow_matrix = new double[max_hmer+1][key.length];
        for (int i = 0 ; i < max_hmer+1; i++) {
            for (int j = 0 ; j < key.length; j++ ){
                flow_matrix[i][j] = fbargs.filling_value;;
            }
        }

        // access qual, convert to ultima representation
        byte[]      quals = samRecord.getBaseQualities();
        byte[]      ti = samRecord.getByteArrayAttribute("ti");
        if ( isReverseStrand() ) {
            reverse(quals, quals.length);
            reverse(ti, ti.length);
        }
        double[]    probs = new double[quals.length];
        for ( int i = 0 ; i < quals.length ; i++ ) {
            double q = quals[i];
            double p = Math.pow(10, -q/10);
            double ultima_p = p*2;

            probs[i] = ultima_p;
        }

        // apply key and qual/ti to matrix
        int     qualOfs = 0;
        for ( int i = 0 ; i < key.length ; i++ ) {
            final byte        run = key[i];
            if ( run <= max_hmer ) {
                flow_matrix[run][i] = (run > 0) ? (1 - probs[qualOfs]) : 1;
                flow_matrix[run][i] = Math.max(MINIMAL_CALL_PROB, flow_matrix[run][i]);
            }
            if ( run != 0 ) {
                if ( quals[qualOfs] != 40 ) {
                    final int     run1 = (ti[qualOfs] == 0) ? (run - 1) : (run + 1);
                    if ( run1 <= max_hmer ) {
                        flow_matrix[run1][i] = probs[qualOfs] / flow_matrix[run][i];
                    }
                    if (run <= max_hmer) {
                        flow_matrix[run][i] /= flow_matrix[run][i]; // for comparison to the flow space - probabilities are normalized by the key's probability
                    }
                }
                qualOfs += run;
            }

        }
    }

    private void readFlowMatrix(String _flow_order) {
        key = getAttributeAsByteArray("kr");
        getKey2Base();

        flow_order = getFlow2Base(_flow_order, key.length);

        flow_matrix = new double[max_hmer+1][key.length];
        for (int i = 0 ; i < max_hmer+1; i++) {
            for (int j = 0 ; j < key.length; j++ ){
                flow_matrix[i][j] = fbargs.filling_value;
            }
        }

        byte [] kh = getAttributeAsByteArray( "kh" );
        int [] kf = getAttributeAsIntArray("kf");
        byte [] kd = getAttributeAsByteArray( "kd");

        byte [] key_kh = key;
        int [] key_kf = new int[key.length];
        for ( int i = 0 ; i < key_kf.length ; i++)
            key_kf[i] = i;
        byte [] key_kd = new byte[key.length];

        kh = ArrayUtils.addAll(kh, key_kh);
        kf = ArrayUtils.addAll(kf, key_kf);
        kd = ArrayUtils.addAll(kd, key_kd);

        quantizeProbs(kd);

        double [] kd_probs = phredToProb(kd);
        if (fbargs.disallow_larger_probs) {
            removeLargeProbs(kd_probs);
        }

        double [] old_kd_probs = ArrayUtils.clone(kd_probs);;

        if (fbargs.remove_longer_than_one_indels) {
            removeLongIndels( key_kh, kh, kf, kd_probs );
        }

        if (fbargs.remove_one_to_zero_probs) {
            removeOneToZeroProbs(key_kh, kh, kf, kd_probs);
        }

        fillFlowMatrix( kh, kf, kd_probs);

        if ((fbargs.lump_probs)) {
            lumpProbs();
        }
        clipProbs();

        if (fbargs.symmetric_indels) {
            smoothIndels(key_kh);
        }
        if (fbargs.only_ins_or_del) {
            reportInsOrDel(key_kh);
        }

        validateSequence();
    }



    public FlowBasedRead(FlowBasedRead other) {
        super(other.samRecord);
        forward_sequence = other.forward_sequence.clone();
        key = other.key.clone();
        flow2base = other.flow2base.clone();
        flow_order = other.flow_order.clone();
        flow_matrix = other.flow_matrix.clone();
        valid_key = other.valid_key;
        direction = other.direction;
        max_hmer = other.max_hmer;
        fbargs = other.fbargs;
    }

    public FlowBasedRead(GATKRead read, String flow_order, int _max_hmer, FlowBasedAlignmentArgumentCollection fbargs) {
        this(read.convertToSAMRecord(null), flow_order, _max_hmer, fbargs);
        this.read = read;

    }

    public String getFlowOrder() {
        return new String(Arrays.copyOfRange(flow_order, 0, Math.min(4,flow_order.length)));
    }

    public int getMaxHmer() {
        return max_hmer;
    }

    public int getNFlows() {
        return key.length;
    }
    public Direction getDirection(){
        return direction;
    }

    private double[] phredToProb(byte [] kq) {
        double [] result = new double[kq.length];
        for (int i = 0 ; i < kq.length; i++ ) {
            result[i] = Math.pow(10, ((double)-kq[i])/10);
        }
        return result;
    }


    private byte[] getForwardSequence(){
        if (!isReverseStrand()) {
            return samRecord.getReadBases();
        } else {
            byte[] result = new byte[samRecord.getReadBases().length];
            System.arraycopy(samRecord.getReadBases(), 0, result, 0, result.length);
            SequenceUtil.reverseComplement(result);
            return result;
        }
    }


    private void fillFlowMatrix(byte [] kh, int [] kf,
                                double [] kd_probs ) {
        for ( int i = 0 ; i < kh.length; i++ ) {
            if (( kh[i] & 0xff )> max_hmer) {
                continue;
            }
            flow_matrix[kh[i] & 0xff][kf[i]] = kd_probs[i];
        }

    }

    private int[] getAttributeAsIntArray(String attributeName) {
        ReadUtils.assertAttributeNameIsLegal(attributeName);
        Object attributeValue = this.samRecord.getAttribute(attributeName);
        if (attributeValue == null) {
            return null;
        } else if (attributeValue instanceof byte[]) {
            byte[] tmp = (byte[]) attributeValue;
            int[] ret = new int[tmp.length];
            for (int i = 0; i < ret.length; i++)
                ret[i] = tmp[i] & 0xff; //converting signed byte to unsigned
            return Arrays.copyOf(ret, ret.length);
        } else if ((attributeValue instanceof int[])) {
            int[] ret = (int[]) attributeValue;
            return Arrays.copyOf(ret, ret.length);
        } else if  (attributeValue instanceof short[]) {
            short [] tmp = (short[]) attributeValue;
            int[] ret = new int[tmp.length];
            for (int i = 0 ; i < tmp.length; i++ )
                ret[i] = tmp[i];
            return Arrays.copyOf(ret, ret.length);
        }else {
            throw new GATKException.ReadAttributeTypeMismatch(attributeName, "integer array");
        }
    }



    private void validateSequence(){
        for (byte b : key) {
            if (b > max_hmer - 1) {
                valid_key = false;
            }
        }
        valid_key = true;
    }

    public boolean is_valid() {
        return valid_key;
    }

    public double getProb(int flow, int hmer) {
        return flow_matrix[Math.min(hmer, max_hmer)][flow];
    }

    public void apply_alignment(){

        if ((getDirection() == Direction.SYNTHESIS) && ( isReverseStrand() )) {
            flipMatrix();
            reverse(key, key.length);
            getKey2Base();
            SequenceUtil.reverseComplement(flow_order);

        }

        boolean isBase = isBaseFormat();
        int[] basePair = {0, 0};
        int[] clip_left_pair = !isBase ? find_left_clipping() : basePair;
        int[] clip_right_pair = !isBase ? find_right_clipping() : basePair;
        int clip_left = clip_left_pair[0];
        int left_hmer_clip = clip_left_pair[1];
        int clip_right = clip_right_pair[0];
        int right_hmer_clip = clip_right_pair[1];

        apply_clipping(clip_left, left_hmer_clip, clip_right, right_hmer_clip);

        setDirection(Direction.REFERENCE);

    }

    private boolean isBaseFormat() {
       return samRecord.hasAttribute("ti");
    }

    public void apply_base_clipping(int clip_left_base, int clip_right_base){
        int[] clip_left_pair = find_left_clipping(clip_left_base);
        int[] clip_right_pair = find_right_clipping(clip_right_base);
        int clip_left = clip_left_pair[0];
        int left_hmer_clip = clip_left_pair[1];
        int clip_right = clip_right_pair[0];
        int right_hmer_clip = clip_right_pair[1];
        if (getLength() - clip_left_base - clip_right_base < MINIMAL_READ_LENGTH) {
            trimmed_to_haplotype = true;
            valid_key=false;
            trim_left_base=clip_left_base;
            trim_right_base = clip_right_base;
        } else {
            apply_clipping(clip_left, left_hmer_clip, clip_right, right_hmer_clip);
            trimmed_to_haplotype = true;
            trim_left_base = clip_left_base;
            trim_right_base = clip_right_base;
        }
    }

    private void apply_clipping(int clip_left, int left_hmer_clip, int clip_right, int right_hmer_clip){
        if ((clip_left < 0) || (clip_right < 0)  || (clip_left >= getKeyLength() ) || ( clip_right >= getKeyLength())) {
            throw new GATKException.ShouldNeverReachHereException("Weird read clip calculated");
            //return 1;
        }

        if ((left_hmer_clip < 0) || (right_hmer_clip < 0)  || (left_hmer_clip >= 14 ) || ( right_hmer_clip >= 14)) {
            throw new GATKException.ShouldNeverReachHereException("Weird read clip calculated");
            //return 1;
        }

        int original_length = key.length;

        key[clip_left]-=left_hmer_clip;
        boolean shift_left = true;
        if ( (clip_left >= 0) || ( left_hmer_clip >= 0 )  ) {
            while (key[clip_left] == 0) {
                clip_left += 1 ;
                shift_left = false;
            }
        }
        key[key.length - clip_right-1] -= right_hmer_clip;
        boolean shift_right = true;
        if ( (clip_right >= 0) || ( right_hmer_clip >= 0 )  ) {
            while (key[original_length - 1- clip_right] == 0) {
                clip_right += 1 ;
                shift_right = false;
            }
        }

        key = Arrays.copyOfRange(key, clip_left, original_length - clip_right);
        getKey2Base();
        flow_order = Arrays.copyOfRange(flow_order, clip_left, original_length - clip_right);

        double [][] new_flow_matrix = new double[flow_matrix.length][original_length - clip_left - clip_right] ;
        for ( int i = 0 ; i < new_flow_matrix.length; i++) {
            new_flow_matrix[i] = Arrays.copyOfRange(flow_matrix[i], clip_left, original_length - clip_right);
        }

        flow_matrix = new_flow_matrix;
        if (shift_left) {
            shiftColumnUp(flow_matrix, 0, left_hmer_clip);
        }

        if (shift_right) {
            shiftColumnUp(flow_matrix, flow_matrix[0].length-1, right_hmer_clip);
        }

    }

    private void reverse(int []a, int n)
    {
        int i, k, t;
        for (i = 0; i < n / 2; i++) {
            t = a[i];
            a[i] = a[n - i - 1];
            a[n - i - 1] = t;
        }

    }
    private void reverse(byte []a, int n)
    {
        int i, k;
        byte t;
        for (i = 0; i < n / 2; i++) {
            t = a[i];
            a[i] = a[n - i - 1];
            a[n - i - 1] = t;
        }

    }

    private void reverse(double []a, int n)
    {
        int i, k;
        double t;
        for (i = 0; i < n / 2; i++) {
            t = a[i];
            a[i] = a[n - i - 1];
            a[n - i - 1] = t;
        }

    }


    private void flipMatrix() {
        for ( int i = 0 ; i < flow_matrix.length; i++) reverse(flow_matrix[i], flow_matrix[i].length);
    }

    private void shiftColumnUp(double[][] matrix, int colnum, int shift) {
        for (int i = 0; i < matrix.length - shift; i ++ ) {
            matrix[i][colnum] = matrix[i+shift][colnum];
        }
        for (int i = matrix.length - shift; i < matrix.length; i ++ ) {
            matrix[i][colnum] = 0;
        }

    }

    public void setDirection( Direction dir ) {
        direction = dir;
    }

    private void getKey2Base() {
        flow2base = new int[key.length];
        flow2base[0] = -1;
        for (int i = 1 ; i < flow2base.length; i++) {
            flow2base[i] = flow2base[i-1] + key[i-1];
        }
    }

    private int[] getKey2Base(byte[] key){
        int[] result = new int[key.length];
        result[0] = -1;
        for (int i = 1 ; i < result.length; i++) {
            result[i] = result[i-1] + key[i-1];
        }
        return result;
    }

    private byte[] getFlow2Base(String flow_order, int expected_length) {
        byte[] result = new byte[expected_length] ;
        for ( int i = 0; i < result.length; i++ ) {
            result[i] = (byte)flow_order.charAt(i%flow_order.length());
        }
        return result;
    }

    private int[] find_left_clipping() {
        List<CigarElement> cigar = getCigarElements();
        int[] result = new int[2];
        if (cigar.size() == 0 ) {
            return result;
        }

        CigarElement start = cigar.get(0);
        if (start.getOperator() != CigarOperator.H) {
            return result;
        }

        int bases_clipped = start.getLength();
        return find_left_clipping(bases_clipped);
    }


    private int[] find_left_clipping(int bases_clipped){
        int[] result = new int[2];
        if (bases_clipped==0){
            return result;
        }
        int index =0 ;
        int stop_clip = 0;
        //System.out.println("Based clipped: "+Integer.toString(bases_clipped));
        for (int i = 0 ; i < flow2base.length; i++ ) {

            if (flow2base[i] + key[i] >= bases_clipped) {
                stop_clip = i;
                break;
            }
        }
        int hmer_clipped = bases_clipped - flow2base[stop_clip] - 1;
        result[0] = stop_clip;
        result[1] = hmer_clipped;
        return result;
    }

    private int[] find_right_clipping() {
        List<CigarElement> cigar = getCigarElements();
        int[] result = new int[2];
        if (cigar.size() == 0 ) {
            result[0] = 0;
            result[1] = 0;
            return result;
        }

        CigarElement end = cigar.get(cigar.size()-1);
        if (end.getOperator() != CigarOperator.H) {
            result[0] = 0;
            result[1] = 0;
            return result;
        }

        int bases_clipped = end.getLength();

        return find_right_clipping(bases_clipped);
    }

    private int[] find_right_clipping(int bases_clipped) {
        int[] result = new int[2];
        if (bases_clipped==0){
            return result;
        }

        int index =0 ;
        int stop_clip = 0;

        byte[] rkey = new byte[key.length];
        for (int i = 0 ; i < key.length; i++ ){
            rkey[i] = key[key.length-1-i];
        }

        int[] rflow2base = getKey2Base(rkey);
        for (int i = 0 ; i < rflow2base.length; i++ ) {
            if (rflow2base[i] + rkey[i] >= bases_clipped) {
                stop_clip = i;
                break;
            }
        }

        int hmer_clipped = bases_clipped - rflow2base[stop_clip] - 1;
        result[0] = stop_clip;
        result[1] = hmer_clipped;
        return result;
    }


    public void writeKey(FileWriter oos)
            throws IOException {
        for (int i = 0; i < key.length; i++)
            oos.write(Byte.toString(key[i])+"\n");

    }

    public void writeMatrix(OutputStreamWriter oos)
            throws IOException {
        DecimalFormat formatter = new DecimalFormat("0.0000", DecimalFormatSymbols.getInstance(Locale.ENGLISH));

        byte[]      bases = samRecord.getReadBases();
        int         basesOfs = 0;
        byte[]      quals = samRecord.getBaseQualities();
        byte[]      ti = samRecord.hasAttribute("ti") ? samRecord.getByteArrayAttribute("ti") : (new byte[key.length*3]);
        if ( isReverseStrand() )
        {
            reverse(quals, quals.length);
            reverse(ti, ti.length);
        }

        for ( int col = 0 ; col < key.length ; col++ ) {
            oos.write("C,R,F,B,Bi,Q,ti\n");
            byte base = (key[col] != 0) ? (basesOfs < bases.length ? bases[basesOfs] : (byte)'?') : (byte)'.';
            String bi = (key[col] != 0) ? Integer.toString(basesOfs) : ".";
            String q = (key[col] != 0) ? Integer.toString(quals[basesOfs]) : ".";
            String Ti = (key[col] != 0) ? Integer.toString(ti[basesOfs]) : ".";
            for (int row = 0; row < flow_matrix.length; row++) {
                String s = formatter.format(flow_matrix[row][col]);
                oos.write(String.format("%d,%d,%d,%c,%s,%s,%s,%s %s\n", col, row, key[col], base, bi, q, Ti, isReverseStrand() ? "r" : ".", s));
            }
            if ( key[col] != 0 )
                basesOfs +=  key[col];
            oos.write("\n");
        }

    }

    public void logMatrix(Logger logger, String msg) {

        ByteArrayOutputStream   os = new ByteArrayOutputStream();
        OutputStreamWriter      osw = new OutputStreamWriter(os);

        logger.debug("logMatrix: " + msg + ", " + getName());

        try {
            writeMatrix(osw);
            for ( String line : (new String(os.toByteArray())).split("\n") )
                logger.debug(line);
            osw.close();;
        } catch (IOException e) {
            e.printStackTrace();;
        }
    }

    public byte [] getFlowOrderArray() {
        return flow_order;
    }

    public int getKeyLength() {
        return key.length;
    }

    public byte[] getKey() {
        return key;
    }

    public int totalKeyBases()  {
        int sum = 0 ;
        for (int i = 0 ; i < key.length; i++){
            sum += key[i];
        }
        return sum;
    }

    public int seqLength(){
        return forward_sequence.length;
    }
    public boolean isTrimmed_to_haplotype() {
        return trimmed_to_haplotype;
    }

    public int getTrimmedStart() {
        return trim_left_base + getStart();
    }
    public int getTrimmedEnd() {
        return getEnd() - trim_right_base;
    }

    //functions that take care of simulating base format
    private void clipProbs() {
        for ( int i = 0 ; i < getMaxHmer(); i++ ) {
            for ( int j =0; j < getNFlows(); j++) {
                if ((flow_matrix[i][j] < fbargs.probability_ratio_threshold) &&
                        (key[j]!=i)) {
                    flow_matrix[i][j] = fbargs.filling_value;
                }
            }
        }
    }

    private void removeLargeProbs(double [] kd_probs) {
        for (int i = 0 ; i < kd_probs.length; i++) {
            if (kd_probs[i] > 1) {
                kd_probs[i] = 1;
            }
        }
    }

    private void removeLongIndels(  byte [] key_kh, byte [] kh, int [] kf, double [] kd_probs ){
        for ( int i = 0 ; i < kd_probs.length; i++ ) {
            if (Math.abs(key_kh[kf[i]] - kh[i]) > 1 ){
                kd_probs[i] = fbargs.filling_value;
            }
        }
    }

    private void removeOneToZeroProbs( byte [] key_kh, byte [] kh, int[] kf, double [] kd_probs) {
        for ( int i = 0 ; i < kd_probs.length; i++ ) {
            if (key_kh[kf[i]]==0 && kh[i]>0 ){
                kd_probs[i] = fbargs.filling_value;
            }
        }
    }




    private void quantizeProbs( byte [] kd_probs ) {
        int nQuants = fbargs.probability_quantization;
        double bin_size = 60/(float)nQuants;
        for ( int i = 0 ; i < kd_probs.length; i++) {
            if (kd_probs[i] <=0)
                continue;
            else {
                kd_probs[i] = (byte)(bin_size * (int)(kd_probs[i]/bin_size)+1);
            }
        }
    }

    private void smoothIndels( byte [] kr ) {
        for ( int i = 0 ; i < kr.length; i++ ){
            byte idx = kr[i];
            if (( idx > 1 ) && ( idx < max_hmer) ) {
                double tmp = (flow_matrix[idx - 1][i] + flow_matrix[idx + 1][i]) / 2;
                flow_matrix[idx - 1][i] = tmp;
                flow_matrix[idx + 1][i] = tmp;
            }
        }
    }

    private void reportInsOrDel( byte [] kr ) {
        for ( int i = 0 ; i < kr.length; i++ ){
            byte idx = kr[i];
            if (( idx > 1 ) && ( idx < max_hmer) ) {
                if ((flow_matrix[idx-1][i] > fbargs.filling_value) && (flow_matrix[idx+1][i] > fbargs.filling_value)) {
                    int fix_cell = flow_matrix[idx-1][i] > flow_matrix[idx+1][i] ? idx+1 : idx-1;
                    int other_cell = flow_matrix[idx-1][i] > flow_matrix[idx+1][i] ? idx-1 : idx+1;
//                    if (fbargs.lump_probs) {
//                       flow_matrix[other_cell][i]+= flow_matrix[fix_cell][i];
//                    }
                    flow_matrix[fix_cell][i] = fbargs.filling_value;
                }
            }
        }
    }

    private void lumpProbs() {

        for (int i = 0; i < getMaxHmer(); i++) {
            for (int j = 0 ; j < getNFlows(); j ++ ) {
                int fkey = key[j];
                if (flow_matrix[i][j]<=fbargs.filling_value) {
                    continue;
                } else {
                    if ( (i - fkey) < -1 ){
                        flow_matrix[fkey-1][j]+=flow_matrix[i][j];
                        flow_matrix[i][j] = fbargs.filling_value;
                    } else if ((i-fkey) > 1) {
                        flow_matrix[fkey+1][j]+=flow_matrix[i][j];
                        flow_matrix[i][j] = fbargs.filling_value;
                    }

                }

            }
        }

    }

}

