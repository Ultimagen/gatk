package org.ultimagenomics.flow_based_read.read;

import org.ultimagenomics.flow_based_read.utils.Direction;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
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
    private final int max_hmer = 9;
    private byte[] flow_order;
    private double[][] flow_matrix;
    private boolean valid_key;
    private Direction direction = Direction.SYNTHESIS;

    public FlowBasedRead(SAMRecord samRecord) {
        super(samRecord);
        this.samRecord = samRecord;
        forward_sequence = getForwardSequence();
        key = getAttributeAsByteArray("ks");
        getKey2Base();
        flow_order = getFlow2Base(getAttributeAsString("KS"), key.length);
        flow_matrix = new double[max_hmer+1][key.length];
        byte [] kq = getAttributeAsByteArray("kq");
        double [] gtr_probs = phredToProb(kq);

        byte [] kh = getAttributeAsByteArray( "kh" );
        int [] kf = getAttributeAsIntArray("kf");
        byte [] kd = getAttributeAsByteArray( "kd");

        double [] kd_probs = phredToProb(kd);
        for ( int i = 0 ; i < kd_probs.length; i++ ) {
            kd_probs[i] = 1 - kd_probs[i];
        }

        fillFlowMatrix( kh, kf, kd_probs, gtr_probs);

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
    }

    public FlowBasedRead(GATKRead read) {
        this(read.convertToSAMRecord(null));
        this.read = read;
    }

    public Direction getDirection(){
        return direction;
    }
    private double[] phredToProb(byte [] kq) {
        double [] result = new double[kq.length];
        for (int i = 0 ; i < kq.length; i++ ) {
            result[i] = 1-Math.pow(10, -(double)kq[i]/10);
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
                                double [] kd_probs, double [] key_probs ) {
        for ( int i = 0 ; i < kh.length; i++ ) {
            flow_matrix[kh[i]&0xff][kf[i]] = kd_probs[i];
        }

        for (int i = 0 ; i < key.length; i++ ) {
             flow_matrix[key[i]][i] = key_probs[i];
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
        return flow_matrix[hmer][flow];
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

        if ((clip_left < 0) || (clip_right < 0)  || (clip_left >= getKeyLength() ) || ( clip_right >= getKeyLength())) {
            throw new GATKException.ShouldNeverReachHereException("Weird read clip calculated");
            //return 1;
        }

        if ((left_hmer_clip < 0) || (right_hmer_clip < 0)  || (left_hmer_clip >= 10 ) || ( right_hmer_clip >= 10)) {
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
        flow2base = Arrays.copyOfRange(flow2base, clip_left, original_length - clip_right);
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
        setDirection(Direction.REFERENCE);

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
}

