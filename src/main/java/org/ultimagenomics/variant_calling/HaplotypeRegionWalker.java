package org.ultimagenomics.variant_calling;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.nio.channels.SeekableByteChannel;
import java.nio.file.Path;
import java.util.LinkedList;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;
import java.util.function.Function;

public class HaplotypeRegionWalker {

    private static final Logger logger = LogManager.getLogger(HaplotypeRegionWalker.class);

    private SamReader       samReader;
    private List<Haplotype> walkerHaplotypes = new LinkedList<>();
    private boolean         reusePreviousResults = false;

    HaplotypeRegionWalker(HaplotypeBasedVariantRecallerArgumentCollection vrArgs) {
        Path samPath = IOUtils.getPath(vrArgs.HAPLOTYPES_BAM_FILE);

        Function<SeekableByteChannel, SeekableByteChannel> cloudWrapper = BucketUtils.getPrefetchingWrapper(40);
        Function<SeekableByteChannel, SeekableByteChannel> cloudIndexWrapper = BucketUtils.getPrefetchingWrapper(40);
        samReader = SamReaderFactory.makeDefault().referenceSequence(vrArgs.REFERENCE_FASTA).open(samPath, cloudWrapper, cloudIndexWrapper);

    }

    void forBest(Locatable queryLoc, Consumer<List<Haplotype>> action) {
        Objects.requireNonNull(action);

        final List<Haplotype> best = new LinkedList<>();

        forEach(queryLoc, haplotypes -> {
            if ( best.size() == 0 ||
                    (fittnessScore(queryLoc, haplotypes) > fittnessScore(queryLoc, best)) ) {
                best.clear();
                best.addAll(haplotypes);
            }
        });

        if ( best.size() != 0 )
            action.accept(best);
    }

    private double fittnessScore(Locatable loc, List<Haplotype> haplotypes) {
        Objects.requireNonNull(haplotypes);
        if ( haplotypes.size() == 0 )
            return 0;
        Locatable   hloc = haplotypes.get(0).getGenomeLocation();

        // determine spacing before and end of loc
        int         before = Math.max(1, loc.getStart() - hloc.getStart());
        int         after = Math.max(1, hloc.getEnd() - loc.getEnd());

        // score reflects closeness to being in the center
        double       score = 1.0 - 2 * Math.abs(0.5 - (double)before / (before + after));

        if ( logger.isDebugEnabled() )
            logger.debug(String.format("loc %s, hloc: %s, before: %d, after: %d, score: %f",
                                                            loc, hloc, before, after, score));

        return score;
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
