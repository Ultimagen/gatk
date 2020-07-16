package org.ultimagenomics.variant_calling;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

import java.util.LinkedList;
import java.util.List;
import java.util.Objects;
import java.util.function.Consumer;

public class HaplotypeRegionWalker {

    private SamReader       samReader;
    private List<Haplotype> haplotypes = new LinkedList<>();

    HaplotypeRegionWalker(HaplotypeBasedVariantRecallerArgumentCollection vrArgs) {

        samReader = SamReaderFactory.makeDefault().referenceSequence(vrArgs.REFERENCE_FASTA).open(vrArgs.HAPLOTYPES_BAM_FILE);
    }

    void forEach(Consumer<List<Haplotype>> action) {
        Objects.requireNonNull(action);


        samReader.forEach(
                record -> {
                    if (record.getReadName().startsWith("HC_") ) {
                        Locatable       loc = new SimpleInterval(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd());
                        if ( (haplotypes.size() != 0) && !haplotypes.get(0).getGenomeLocation().equals(loc) ) {
                            action.accept(haplotypes);
                            haplotypes.clear();
                        }
                        Haplotype       haplotype = new Haplotype(record.getReadBases(), loc);
                        haplotype.setCigar(record.getCigar());
                        haplotypes.add(haplotype);
                    }
                }
        );
        if ( haplotypes.size() != 0 )
            action.accept(haplotypes);
    }

}
