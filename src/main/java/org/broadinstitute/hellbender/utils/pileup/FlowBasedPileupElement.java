package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.samtools.CigarElement;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.GATKRead;

public class FlowBasedPileupElement {
    private FlowBasedRead read;
    private PileupElement pe;
    public FlowBasedPileupElement(final FlowBasedRead read,
                                  final PileupElement pe) {
        this.read = read;
        this.pe = pe;
    }

    public final int getFlowNum() {
        return read.getBase2Flow(pe.getOffset());
    }

    public final byte getFlowNuc() {
        return read.getNucForFlow(read.getBase2Flow(pe.getOffset()));
    }

    public final int getHpolCall() {
        return read.getKey()[read.getBase2Flow(pe.getOffset())];
    }
    public final double[] getCallProbs(){
        return read.getHmerProbs(read.getBase2Flow(pe.getOffset()));
    }
}
