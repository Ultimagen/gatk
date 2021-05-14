package org.ultimagen.flowBasedRead.alignment;

import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

public class AlignmentThreadingUtils {

    static private ConcurrentMap<Long, SmithWatermanAligner> threadAligners = new ConcurrentHashMap<>();
    static private boolean singleAligner = true;

    public static SmithWatermanAligner getSimilarAligner(SmithWatermanAligner aligner) {

        SmithWatermanAligner.Implementation implementation = aligner.getClass().getName().contains("Intel")
                ? SmithWatermanAligner.Implementation.FASTEST_AVAILABLE
                : SmithWatermanAligner.Implementation.JAVA;
        return SmithWatermanAligner.getAligner(implementation);
    }

    public static SmithWatermanAligner getSimilarAlignerForCurrentThread(SmithWatermanAligner aligner) {

        if ( singleAligner || Thread.currentThread().getName().equals("main") )
            return aligner;
        long                        id = Thread.currentThread().getId();
        SmithWatermanAligner        threadAligner = threadAligners.get(id);
        if ( threadAligner == null ) {
            threadAligners.put(id, threadAligner = getSimilarAligner(aligner));
        }
        return threadAligner;
    }
}
