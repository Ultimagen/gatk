package org.ultimagenomics.flow_based_read.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

import java.util.ArrayList;

public class FlowBasedHaplotype  extends Allele {
    private static final long serialVersionUID = 42L;
    private static int N_ASCII=78;
    private byte [] key;
    private byte [] r_key;
    private int[] flow2base;
    private int[] r_flow2base;
    private String flow_order;
    private Locatable genomeLoc;
    private Cigar cigar;
    private int ref_end;
    private byte[] flow_order_array;

    public int getKeyLength() {
        return key.length;
    }

    public byte[] getKey() {
        return key;
    }


    public FlowBasedHaplotype(Haplotype source_haplotype, String flow_order){
        this(source_haplotype, flow_order, 0 );
    }

    public FlowBasedHaplotype(Haplotype source_haplotype, String flow_order, int clipping){
        super(source_haplotype.getBases(), source_haplotype.isReference());
        key = base2key(source_haplotype.getBases(), flow_order, 1000 );
        this.flow_order = flow_order;
        genomeLoc = source_haplotype.getGenomeLocation();
        cigar = source_haplotype.getCigar();
        flow2base = getKey2Base(key);
        r_key = key.clone();
        reverse(r_key, r_key.length);
        r_flow2base = getKey2Base(r_key);
        flow_order_array = getFlow2Base(flow_order, key.length);
    }

    private byte[] getFlow2Base(String flow_order, int expected_length) {
        byte[] result = new byte[expected_length] ;
        for ( int i = 0; i < result.length; i++ ) {
            result[i] = (byte)flow_order.charAt(i%flow_order.length());
        }
        return result;
    }



    private int[] getKey2Base(byte[] key){
        int[] result = new int[key.length];
        result[0] = -1;
        for (int i = 1 ; i < result.length; i++) {
            result[i] = result[i-1] + key[i-1];
        }
        return result;
    }

    static private byte[] base2key(byte[] bases, String flow_order, int clipping){
        ArrayList<Byte> result = new ArrayList<>();
        byte[] flow_order_bytes = flow_order.getBytes();
        int loc = 0;
        int flow_number = 0 ;
        int period = flow_order_bytes.length;
        while ( loc < bases.length ) {
            byte flowBase = flow_order_bytes[flow_number%period];
            if ((bases[loc]!=flowBase) && ( bases[loc]!=N_ASCII)) {
                result.add((byte) 0);
            } else {
                int count = 0;
                while ( ( loc < bases.length) && ((bases[loc]==flowBase) || (bases[loc]==N_ASCII)) ){
                    loc++;
                    count ++;
                }
                result.add(count<clipping ? (byte)count : (byte)clipping);
            }
            flow_number++;
        }
        byte[] ret = new byte[result.size()];
        for (int i = 0; i < result.size(); i++) {
            ret[i] = result.get(i);
        }
        return ret;
    }

    public int getStart(){
        return genomeLoc.getStart();
    }

    public int getEnd() {
        return genomeLoc.getEnd();
    }

    public String getChr() {
        return genomeLoc.getContig();
    }

    public Cigar getCigar(){
        return cigar;
    }

    public int[] find_left_clipping(int base_clipping) {
        int [] result = new int[2];
        if (base_clipping == 0 ){
            return result;
        }


        int index =0 ;
        int stop_clip = 0;
        //System.out.println("Based clipped: "+Integer.toString(bases_clipped));
        for (int i = 0 ; i < flow2base.length; i++ ) {

            if (flow2base[i] + key[i] >= base_clipping) {
                stop_clip = i;
                break;
            }
        }
        int hmer_clipped = base_clipping - flow2base[stop_clip] - 1;
        result[0] = stop_clip;
        result[1] = hmer_clipped;
        return result;
    }

    public int[] find_right_clipping(int base_clipping) {
        int [] result = new int[2];
        if (base_clipping == 0 ){
            return result;
        }


        int index = 0;
        int stop_clip = 0;

        for (int i = 0 ; i < r_flow2base.length; i++ ) {
            if (r_flow2base[i] + r_key[i] >= base_clipping) {
                stop_clip = i;
                break;
            }
        }

        int hmer_clipped = base_clipping - r_flow2base[stop_clip] - 1;
        result[0] = stop_clip;
        result[1] = hmer_clipped;
        return result;
    }

    private void reverse(int [] a, int n)
    {
        int i, k, t;
        for (i = 0; i < n / 2; i++) {
            t = a[i];
            a[i] = a[n - i - 1];
            a[n - i - 1] = t;
        }

    }
    private void reverse(byte [] a, int n)
    {
        int i, k;
        byte t;
        for (i = 0; i < n / 2; i++) {
            t = a[i];
            a[i] = a[n - i - 1];
            a[n - i - 1] = t;
        }

    }

    private void reverse(double [] a, int n)
    {
        int i, k;
        double t;
        for (i = 0; i < n / 2; i++) {
            t = a[i];
            a[i] = a[n - i - 1];
            a[n - i - 1] = t;
        }

    }

    public byte [] getFlowOrderArray() {
        return flow_order_array;
    }


}
