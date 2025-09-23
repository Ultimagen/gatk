package org.broadinstitute.hellbender.utils.pileup;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.stream.Collectors;

public class FlowBasedReadPileup extends ReadPileup {
    //private List<FlowBasedPileupElement> fbPileupElements;
    private final List<FlowBasedRead> fbReads;
    private List<FlowBasedPileupElement> fbPile;
    public FlowBasedReadPileup(Locatable loc, List<FlowBasedRead> reads)
    {
        super(loc, reads.stream().map(r -> (GATKRead) r).collect(Collectors.toList()));
        fbReads = reads;
        fbPile = this.pileupElements.stream()
                .map(pe -> new FlowBasedPileupElement((FlowBasedRead)pe.getRead(), pe))
                .collect(Collectors.toList());
    }
}
