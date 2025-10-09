package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.FlowBasedReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class FlowBasedAlignmentLikelihoodEngineTestUtils {

    public static AlleleLikelihoods<GATKRead, Haplotype> computeReadLikelihoods(
            final List<Haplotype> haplotypeList, final List<GATKRead> reads,
            final boolean filterPoorly, final SAMFileHeader hdr,
            final FlowBasedAlignmentLikelihoodEngine engine) {


        final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(haplotypeList);
        final ArrayList<String> _sampList = new ArrayList<>();
        _sampList.add("HG001");
        final SampleList samples = new IndexedSampleList(_sampList);

        // Add likelihoods for each sample's reads to our result
        final HashMap<String, List<GATKRead>> perSampleReadList = new HashMap<>();
        final ArrayList<GATKRead> freads = new ArrayList<>();

        final FlowBasedReadUtils.ReadGroupInfo rgInfo = FlowBasedReadUtils.getReadGroupInfo(hdr, reads.get(0));
        final String flowOrder = rgInfo.flowOrder.substring(0, 4);
        final int maxHmer = rgInfo.maxClass;

        for ( final GATKRead r : reads ) {
            freads.add(new FlowBasedRead(r, flowOrder, maxHmer, new FlowBasedAlignmentArgumentCollection()));
        }
        perSampleReadList.put("HG001", freads);

        final AlleleLikelihoods<GATKRead, Haplotype> result = new AlleleLikelihoods<>(samples, haplotypes, perSampleReadList);
        final int sampleCount = result.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            engine.computeReadLikelihoods(result.sampleMatrix(i), flowOrder);
        }

        result.normalizeLikelihoods(engine.getLog10globalReadMismappingRate(), engine.isSymmetricallyNormalizeAllelesToReference());
        if ( filterPoorly ) {
            result.filterPoorlyModeledEvidence(engine.log10MinTrueLikelihood(engine.getExpectedErrorRatePerBase(), false));
        }

        return result;
    }
}
