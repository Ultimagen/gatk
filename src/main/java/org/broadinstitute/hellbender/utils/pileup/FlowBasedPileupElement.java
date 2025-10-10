package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.samtools.CigarElement;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class FlowBasedPileupElement {
    private FlowBasedRead read;
    private PileupElement pe;
    private boolean offsetDeletions;
    public FlowBasedPileupElement(final FlowBasedRead read,
                                  final PileupElement pe,
                                  final boolean offsetDeletions) {
        this.read = read;
        this.pe = pe;
        this.offsetDeletions = offsetDeletions;
    }
    private final int getBaseOffsetWithDeletion(PileupElement pe){
        if (pe.isDeletion() && offsetDeletions) {
            return pe.getOffset()+1;
        } else {
            return pe.getOffset();
        }
    }
    public final int getFlowNum() {
        return read.getBase2Flow(getBaseOffsetWithDeletion(pe));
    }

    public final byte getFlowNuc() {
        return read.getNucForFlow(read.getBase2Flow(getBaseOffsetWithDeletion(pe)));
    }

    public final int getHpolCall() {
        return read.getKey()[read.getBase2Flow(getBaseOffsetWithDeletion(pe))];
    }
    public final double[] getCallProbs(){
        return read.getHmerProbs(read.getBase2Flow(getBaseOffsetWithDeletion(pe)));
    }
}
