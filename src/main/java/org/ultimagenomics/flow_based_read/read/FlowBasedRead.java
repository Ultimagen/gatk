package org.ultimagenomics.flow_based_read.read;

import htsjdk.samtools.*;
import org.apache.commons.lang.ArrayUtils;
import org.ultimagenomics.flow_based_read.utils.Direction;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

public class FlowBasedRead extends SAMRecordToGATKReadAdapter implements GATKRead, FlowBasedReadInterface, Serializable {

    private static final long serialVersionUID = 42L;
    private GATKRead read = null;
    private SAMRecord samRecord;
    private byte[] forward_sequence;
    private byte[] key;
    private int [] flow2base;
    private int max_hmer = 13;
    private byte[] flow_order;
    private double[][] flow_matrix;
    private boolean valid_key;
    private Direction direction = Direction.SYNTHESIS;
    private boolean trimmed_to_haplotype = false;
    private int trim_left_base = 0 ;
    private int trim_right_base = 0 ;
    private final int MINIMAL_READ_LENGTH = 10; // check if this is the right number
    private final double ERROR_PROB=1e-4;

    public FlowBasedRead(SAMRecord samRecord, int _max_hmer) {
        super(samRecord);
        max_hmer = _max_hmer;
        this.samRecord = samRecord;

        forward_sequence = getForwardSequence();
        key = getAttributeAsByteArray("kr");
        getKey2Base();

        flow_order = getFlow2Base(getAttributeAsString("KS"), key.length);
        flow_matrix = new double[max_hmer+1][key.length];

        // we fill all values in the matrix with non-zeros to not discard any read
        for (int i = 0 ; i < max_hmer+1; i++) {
            for (int j = 0 ; j < key.length; j++ ){
                flow_matrix[i][j] = ERROR_PROB;
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

        double [] kd_probs = phredToProb(kd);

        fillFlowMatrix( kh, kf, kd_probs);

        validateSequence();
    }

    public FlowBasedRead(SAMRecord samRecord, String _flow_order, int _max_hmer) {
        super(samRecord);
        max_hmer = _max_hmer;
        this.samRecord = samRecord;
        forward_sequence = getForwardSequence();
        key = getAttributeAsByteArray("kr");
        getKey2Base();

        flow_order = getFlow2Base(_flow_order, key.length);

        flow_matrix = new double[max_hmer+1][key.length];
        for (int i = 0 ; i < max_hmer+1; i++) {
            for (int j = 0 ; j < key.length; j++ ){
                flow_matrix[i][j] = ERROR_PROB;
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

        double [] kd_probs = phredToProb(kd);

        fillFlowMatrix( kh, kf, kd_probs);

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
    }

    public FlowBasedRead(GATKRead read, String flow_order, int _max_hmer) {
        this(read.convertToSAMRecord(null), flow_order, _max_hmer);
        this.read = read;
    }

    public String getFlowOrder() {
        return new String(Arrays.copyOfRange(flow_order, 0, Math.min(4,flow_order.length)));
    }

    public int getMaxHmer() {
        return max_hmer;
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

        int[] clip_left_pair = find_left_clipping();
        int[] clip_right_pair = find_right_clipping();
        int clip_left = clip_left_pair[0];
        int left_hmer_clip = clip_left_pair[1];
        int clip_right = clip_right_pair[0];
        int right_hmer_clip = clip_right_pair[1];

        apply_clipping(clip_left, left_hmer_clip, clip_right, right_hmer_clip);

        setDirection(Direction.REFERENCE);

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

    public void writeMatrix(FileWriter oos)
            throws IOException {
        for (double[] flowMatrix : flow_matrix)
            for (int j = 0; j < flowMatrix.length; j++) {
                DecimalFormat formatter = new DecimalFormat("0.0000", DecimalFormatSymbols.getInstance(Locale.ENGLISH));
                String s = formatter.format(flowMatrix[j]);
                oos.write(s);
                oos.write("\n");
            }
    }

    public byte [] getFlowOrderArray() {
        return flow_order;
    }

    public int getKeyLength() {
        return key.length;
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
}

