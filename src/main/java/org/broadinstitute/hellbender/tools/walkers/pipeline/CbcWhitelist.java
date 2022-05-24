package org.broadinstitute.hellbender.tools.walkers.pipeline;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;

public class CbcWhitelist {

    private static final Logger logger = LogManager.getLogger(CbcWhitelist.class);

    enum ResultType {
        NOT_FOUND,
        EXACT,
        DEL_CORRECTED,
        INS_CORRECTED,
        AMBIGUOUS
    }

    static class Result {
        public final ResultType type;
        public final String correction;

        Result(final ResultType type, final String correction) {
            this.type = type;
            this.correction = correction;
        }

        Result(final ResultType type) {
            this(type, null);
        }

    }

    final Set<String>  seqs = new LinkedHashSet<>();
    final Map<String, String> indelCorrections = new LinkedHashMap<>();
    final Set<String> indelAmbiguous = new LinkedHashSet<>();
    final Set<String> observedSeqs = new LinkedHashSet<>();

    static final Result notFoundResult = new Result(ResultType.NOT_FOUND);
    static final Result exactResult = new Result(ResultType.EXACT);
    static final Result ambiguousResult = new Result(ResultType.AMBIGUOUS);

    public CbcWhitelist(final GATKPath path) {

        // read the file
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(path.getInputStream())) ) {
            String line;
            while ( (line = reader.readLine()) != null ) {
                seqs.add(line);
                if ( (seqs.size() % 1000000) == 0 ) {
                    logger.info("read " + seqs.size() + " whitelist sequences");
                }
            }
        } catch (IOException e) {
            throw new GATKException("failed to open whitelist file", e);
        }
        logger.info("read " + seqs.size() + " whitelist sequences from " + path);
    }

    private Set<String> buildIndelPermutations(final String seq, final boolean permutateInserts, final int resultSizeLimit) {

        // build hmers
        BaseUtils.HmerIterator      iter = new BaseUtils.HmerIterator(seq.getBytes());
        Set<String>                 indels = new LinkedHashSet<>();
        int                         hmerOfs = 0;
        while ( iter.hasNext() && (indels.size() < resultSizeLimit) ) {
            // get hmer
            final Pair<Byte, Integer> hmer = iter.next();
            final int hmerLength = hmer.getRight();

            if ( permutateInserts ) {
                // ins?
                if (hmerLength < FlowBasedRead.MAX_CLASS) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(seq, 0, hmerOfs + hmerLength);
                    sb.append(seq, hmerOfs + hmerLength - 1, hmerOfs + hmerLength); // insert base
                    sb.append(seq, hmerOfs + hmerLength, seq.length());
                    final String indel = sb.toString();
                    if ( seqs.contains(indel) ) {
                        indels.add(sb.toString());
                    }
                }
            } else {
                // del?
                if (hmerLength >= 2) {

                    StringBuilder sb = new StringBuilder();
                    sb.append(seq, 0, hmerOfs + hmerLength - 1); // delete base
                    sb.append(seq, hmerOfs + hmerLength, seq.length());
                    final String indel = sb.toString();
                    if ( seqs.contains(indel) ) {
                        indels.add(sb.toString());
                    }
                }
            }

            // update (next) hmer offset
            hmerOfs += hmerLength;
        }

        return indels;
    }

    private void processObservedSequence(final String seq, final boolean permutateInserts) {

        if ( !observedSeqs.contains(seq) ) {
            Set<String>         permutations = buildIndelPermutations(seq, permutateInserts, 2);
            if ( permutations.size() == 1 ) {
                indelCorrections.put(seq, permutations.iterator().next());
            } else if ( permutations.size() > 1 ) {
                indelAmbiguous.add(seq);
            }
            observedSeqs.add(seq);
        }
    }

    public Result matchOrCorrect(final byte[] bytes, final int start, final int length) {

        // check for an exact match
        if ( isInWhitelist(bytes, start, length) ) {
            return exactResult;
        }

        // check for del
        final String delSeq = new String(Arrays.copyOfRange(bytes, start, start + length - 1));
        processObservedSequence(delSeq, true);
        final String delCorrection;
        if ( indelAmbiguous.contains(delSeq) ) {
            return ambiguousResult;
        } else if (indelCorrections.containsKey(delSeq) ) {
            delCorrection = indelCorrections.get(delSeq);
        } else {
            delCorrection = null;
        }

        // check for ins
        final String insSeq = new String(Arrays.copyOfRange(bytes, start, start + length + 1));
        processObservedSequence(delSeq, false);
        final String insCorrection;
        if ( indelAmbiguous.contains(insSeq) ) {
            return ambiguousResult;
        } else if (indelCorrections.containsKey(insSeq) ) {
            insCorrection = indelCorrections.get(insSeq);
        } else {
            insCorrection = null;
        }

        // build result
        if ( delCorrection != null ) {
            return (insCorrection != null) ? ambiguousResult :  new Result(ResultType.DEL_CORRECTED, delCorrection);
        } else {
            return (insCorrection != null) ? new Result(ResultType.INS_CORRECTED, insCorrection) : notFoundResult;
        }
    }

    public boolean isInWhitelist(final String seq) {
        return seqs.contains(seq);
    }

    public boolean isInWhitelist(final byte[] bytes, final int start, final int length) {
        return isInWhitelist(new String(Arrays.copyOfRange(bytes, start, start + length)));
    }
}
