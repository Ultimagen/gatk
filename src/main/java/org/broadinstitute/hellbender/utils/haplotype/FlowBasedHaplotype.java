package org.broadinstitute.hellbender.utils.haplotype;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.utils.read.FlowBasedKeyCodec;
import org.broadinstitute.hellbender.utils.read.FlowBasedReadUtils;

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
        key = FlowBasedKeyCodec.base2key(sourceHaplotype.getBases(), flowOrder);
        genomeLoc = sourceHaplotype.getGenomeLocation();
        cigar = sourceHaplotype.getCigar();
        flow2base = FlowBasedKeyCodec.getKey2Base(key);
        rKey = key.clone();
        ArrayUtils.reverse(rKey);
        rFlow2Base = FlowBasedKeyCodec.getKey2Base(rKey);
        flowOrderArray = FlowBasedKeyCodec.getFlow2Base(flowOrder, key.length);
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

    public int[] findLeftClipping(final int baseClipping) {
        return FlowBasedReadUtils.findLeftClipping(baseClipping, flow2base, key);
    }

    public int[] findRightClipping(final int baseClipping) {
        return FlowBasedReadUtils.findRightClipping(baseClipping, rFlow2Base, rKey);
    }

    public byte [] getFlowOrderArray() {
        return flowOrderArray;
    }
}
