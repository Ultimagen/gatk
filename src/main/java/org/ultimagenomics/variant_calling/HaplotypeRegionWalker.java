package org.ultimagenomics.variant_calling;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.ultimagenomics.flow_based_read.read.FlowBasedRead;

import java.util.LinkedList;
import java.util.List;
import java.util.Objects;
import java.util.function.Consumer;

public class HaplotypeRegionWalker {

    private SamReader       samReader;
    private List<Haplotype> walkerHaplotypes = new LinkedList<>();

    HaplotypeRegionWalker(HaplotypeBasedVariantRecallerArgumentCollection vrArgs) {

        samReader = SamReaderFactory.makeDefault().referenceSequence(vrArgs.REFERENCE_FASTA).open(vrArgs.HAPLOTYPES_BAM_FILE);
    }

    void forEach(Locatable queryLoc, Consumer<List<Haplotype>> action) {
        Objects.requireNonNull(action);

        // use last results? (note that this can be problematic if a vc is inside two haplotype areas)
        if ( (walkerHaplotypes.size() != 0) && walkerHaplotypes.get(0).getGenomeLocation().contains(queryLoc) ) {
            action.accept(walkerHaplotypes);
        } else {

            // must query
            walkerHaplotypes.clear();
            SAMRecordIterator iter = samReader.query(queryLoc.getContig(), queryLoc.getStart(), queryLoc.getEnd(), false);
            iter.forEachRemaining(record -> {
                if (isHaplotypeRecord(record)) {
                    Locatable loc = new SimpleInterval(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd());
                    if ((walkerHaplotypes.size() != 0) && !walkerHaplotypes.get(0).getGenomeLocation().equals(loc)) {
                        action.accept(walkerHaplotypes);
                        walkerHaplotypes.clear();
                    }
                    walkerHaplotypes.add(buildHaplotype(record));
                }
            });
            if (walkerHaplotypes.size() != 0)
                action.accept(walkerHaplotypes);
            iter.close();
        }
    }

    private boolean isHaplotypeRecord(SAMRecord record) {
        return record.getReadName().startsWith("HC_");
    }

    private Haplotype buildHaplotype(SAMRecord record) {

        Locatable loc = new SimpleInterval(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd());
        Haplotype haplotype = new Haplotype(record.getReadBases(), loc);
        haplotype.setCigar(record.getCigar());
        return haplotype;
    }
}
