package org.broadinstitute.hellbender.tools.walkers.groundtruth;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.Tuple;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

public class AncestralContigLocationTranslator {

    // locals
    final private GATKPath                                    basePath;
    final private Map<String, SingleFileLocatonTranslator>    translators = new LinkedHashMap<>();

    AncestralContigLocationTranslator(GATKPath basePath) {
        this.basePath = basePath;
    }

    Tuple<SimpleInterval, SimpleInterval> translate(final Locatable loc) throws IOException {
        return new Tuple<>(translate(GroundTruthConstants.C_MATERNAL, loc),
                translate(GroundTruthConstants.C_PATERNAL, loc));
    }

    SimpleInterval translate(final String ancestor, final Locatable loc) throws IOException {

        int         start = translate(ancestor, loc.getContig(), loc.getStart());
        int         end = translate(ancestor, loc.getContig(), loc.getEnd());

        if ( end > start ) {
            return new SimpleInterval(loc.getContig() + "_" + ancestor, start, end);
        } else {
            throw new LocationTranslationException("location " + loc + " failed to translate for " + ancestor + ", start:" + start + " ,end:" + end);
        }
    }

    int translate(final String ancestor, final String contig, final int from) throws IOException {

        // check-for/create translator
        final String                          key = ancestor + "." + contig + ".csv";
        if ( !translators.containsKey(key) ) {
            final GATKPath        path = new GATKPath(basePath.getURIString() + key);
            translators.put(key, new SingleFileLocatonTranslator(path));
        }

        // translate
        return translators.get(key).translate(from);
    }
}
