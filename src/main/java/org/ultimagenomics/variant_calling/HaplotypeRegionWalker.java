package org.ultimagenomics.variant_calling;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.nio.channels.SeekableByteChannel;
import java.nio.file.Path;
import java.util.LinkedList;
import java.util.List;
import java.util.Objects;
import java.util.function.Consumer;
import java.util.function.Function;

public class HaplotypeRegionWalker {

    private SamReader       samReader;
    private List<Haplotype> walkerHaplotypes = new LinkedList<>();
    private boolean         reusePreviousResults = false;

    HaplotypeRegionWalker(HaplotypeBasedVariantRecallerArgumentCollection vrArgs) {
        Path samPath = IOUtils.getPath(vrArgs.HAPLOTYPES_BAM_FILE);

        Function<SeekableByteChannel, SeekableByteChannel> cloudWrapper = BucketUtils.getPrefetchingWrapper(40);
        Function<SeekableByteChannel, SeekableByteChannel> cloudIndexWrapper = BucketUtils.getPrefetchingWrapper(40);
        samReader = SamReaderFactory.makeDefault().referenceSequence(vrArgs.REFERENCE_FASTA).open(samPath, cloudWrapper, cloudIndexWrapper);

    }

    void forEach(Locatable queryLoc, Consumer<List<Haplotype>> action) {
        Objects.requireNonNull(action);

        // use last results? (note that this can be problematic if a vc is inside two haplotype areas)
        if ( reusePreviousResults
                && (walkerHaplotypes.size() != 0)
                && walkerHaplotypes.get(0).getGenomeLocation().contains(queryLoc) ) {
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
