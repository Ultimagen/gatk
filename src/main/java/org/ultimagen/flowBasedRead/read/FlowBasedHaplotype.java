package org.ultimagen.flowBasedRead.read;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

/**
 * Haplotype that also keeps information on the flow space @see FlowBasedRead
 * Haplotype can't be extended, so this extends Allele
 */
public class FlowBasedHaplotype  extends Allele {
    private static final long serialVersionUID = 42L;
    private int [] key;
    private int [] rKey;
    private int[] flow2base;
    private int[] rFlow2Base;
    private Locatable genomeLoc;
    private Cigar cigar;
    private byte[] flowOrderArray;

    /* Create flow based haplotype from the haplotype */
    public FlowBasedHaplotype(final Haplotype sourceHaplotype, final String flowOrder){
        super(sourceHaplotype.getBases(), sourceHaplotype.isReference());
        key = FlowBasedRead.base2key(sourceHaplotype.getBases(), flowOrder);
        genomeLoc = sourceHaplotype.getGenomeLocation();
        cigar = sourceHaplotype.getCigar();
        flow2base = FlowBasedRead.getKey2Base(key);
        rKey = key.clone();
        reverse(rKey, rKey.length);
        rFlow2Base = FlowBasedRead.getKey2Base(rKey);
        flowOrderArray = FlowBasedRead.getFlow2Base(flowOrder, key.length);
    }


    public int getKeyLength() {
        return key.length;
    }

    public int[] getKey() {
        return key;
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

    /* clips flows from the left to clip the input number of bases
    Needed to trim the haplotype to the read
    Returns number of flows to remove and the change in the left most remaining flow if necessary
     */
    public int[] findLeftClipping(final int baseClipping) {
        final int [] result = new int[2];
        if (baseClipping == 0 ){
            return result;
        }

        int stopClip = 0;
        for (int i = 0 ; i < flow2base.length; i++ ) {

            if (flow2base[i] + key[i] >= baseClipping) {
                stopClip = i;
                break;
            }
        }
        final int hmerClipped = baseClipping - flow2base[stopClip] - 1;
        result[0] = stopClip;
        result[1] = hmerClipped;
        return result;
    }

    /* clips flows from the right to trim the input number of bases
    Returns number of flows to remove and the change in the right most flow.
     */

    public int[] findRightClipping(final int baseClipping) {
        final int [] result = new int[2];
        if (baseClipping == 0 ){
            return result;
        }


        int stopClip = 0;

        for (int i = 0; i < rFlow2Base.length; i++ ) {
            if (rFlow2Base[i] + rKey[i] >= baseClipping) {
                stopClip = i;
                break;
            }
        }

        final int hmerClipped = baseClipping - rFlow2Base[stopClip] - 1;
        result[0] = stopClip;
        result[1] = hmerClipped;
        return result;
    }

    private void reverse(final int [] a, final int n)
    {
        int i, k, t;
        for (i = 0; i < n / 2; i++) {
            t = a[i];
            a[i] = a[n - i - 1];
            a[n - i - 1] = t;
        }

    }

    public byte [] getFlowOrderArray() {
        return flowOrderArray;
    }

    //check if two haplotypes are equal up to a single hmer change
    public boolean equalUpToHmerChange(final FlowBasedHaplotype other ) {
        if (other.getKeyLength() != getKeyLength() ){
            return false;
        }
        int diffCounts = 0;
        final int [] otherkey = other.getKey();
        for (int i = 0 ; i < getKeyLength(); i ++ ) {
            if (((key[i]==0) && (otherkey[i]!=0)) ||
                    ((key[i]!=0) && ( otherkey[i]==0))) {
                return false;
            }

            if (Math.abs(key[i]-otherkey[i]) > 0 ){
                diffCounts ++;
            }
            if (diffCounts > 1 ) {
                return false;
            }
        }
        return true;
    }
}
