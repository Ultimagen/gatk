package org.ultimagenomics.variant_calling;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Objects;
import java.util.function.Consumer;

public class HaplotypeRegionWalker {

    private SamReader       samReader;
    private Locatable       lastLoc;

    HaplotypeRegionWalker(HaplotypeBasedVariantRecallerArgumentCollection vrArgs) {

        samReader = SamReaderFactory.makeDefault().referenceSequence(vrArgs.REFERENCE_FASTA).open(vrArgs.HAPLOTYPES_BAM_FILE);
        lastLoc = null;
    }

    void forEach(Consumer<Locatable> action) {
        Objects.requireNonNull(action);


        samReader.forEach(
                record -> {
                    if (record.getReadName().startsWith("HC_") ) {
                        Locatable       loc = new SimpleInterval(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd());
                        if ( lastLoc == null || !loc.equals(lastLoc) ) {
                            lastLoc = loc;
                            action.accept(loc);
                        }
                    }
                }
        );
    }

}
