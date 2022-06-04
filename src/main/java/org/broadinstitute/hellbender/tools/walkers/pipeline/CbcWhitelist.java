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
        SNP_CORRECTED,
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

    final boolean supportsSnp;
    final Set<String>  seqs = new LinkedHashSet<>();
    final Map<String, String> indelCorrections = new LinkedHashMap<>();
    final Set<String> indelAmbiguous = new LinkedHashSet<>();
    final Set<String> observedSeqs = new LinkedHashSet<>();

    static final Result notFoundResult = new Result(ResultType.NOT_FOUND);
    static final Result exactResult = new Result(ResultType.EXACT);
    static final Result ambiguousResult = new Result(ResultType.AMBIGUOUS);

    public CbcWhitelist(final GATKPath path, final boolean supportsSnp) {

        this.supportsSnp = supportsSnp;

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

    private Set<String> buildSnpPermutations(final String seq, final int resultSizeLimit) {

        final byte[] allBases = {'A', 'G', 'C', 'T'};
        final Set<String> snps = new LinkedHashSet<>();

        // build hmers
        for ( int ofs = 0 ; ofs < seq.length() ; ofs++ ) {
            final byte base = (byte)seq.charAt(ofs);
            for ( byte snp : allBases ) {
                if ( snp != base ) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(seq, 0, ofs);
                    sb.append((char)snp);
                    sb.append(seq, ofs + 1, seq.length());
                    final String snpSeq = sb.toString();
                    if ( seqs.contains(snpSeq) ) {
                        snps.add(sb.toString());
                        if ( snps.size() >= resultSizeLimit )
                            return snps;
                    }
                }
            }
        }

        return snps;
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

    private void processObservedSequence(final String seq) {

        if ( !observedSeqs.contains(seq) ) {
            Set<String>         permutations = buildSnpPermutations(seq, 2);
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
        final String delCorrection;
        {
            final String seq = new String(Arrays.copyOfRange(bytes, start, start + length - 1));
            processObservedSequence(seq, true);
            if (indelAmbiguous.contains(seq)) {
                return ambiguousResult;
            } else if (indelCorrections.containsKey(seq)) {
                delCorrection = indelCorrections.get(seq);
            } else {
                delCorrection = null;
            }
        }

        // check for ins
        final String insCorrection;
        {
            final String seq = new String(Arrays.copyOfRange(bytes, start, start + length + 1));
            processObservedSequence(seq, false);
            if (indelAmbiguous.contains(seq)) {
                return ambiguousResult;
            } else if (indelCorrections.containsKey(seq)) {
                insCorrection = indelCorrections.get(seq);
            } else {
                insCorrection = null;
            }
        }

        // already ambig?
        if ( delCorrection != null && insCorrection != null ) {
            return ambiguousResult;
        }

        // look for snp
        final String snpCorrection;
        if ( supportsSnp ) {
            final String seq = new String(Arrays.copyOfRange(bytes, start, start + length));
            processObservedSequence(seq);
            if (indelAmbiguous.contains(seq)) {
                return ambiguousResult;
            } else if (indelCorrections.containsKey(seq)) {
                snpCorrection = indelCorrections.get(seq);
            } else {
                snpCorrection = null;
            }
        } else {
            snpCorrection = null;
        }

        // build result
        if ( delCorrection != null ) {
            if ( insCorrection != null ) {
                return ambiguousResult;
            } else {
                return (snpCorrection != null) ? ambiguousResult : new Result(ResultType.DEL_CORRECTED, delCorrection);
            }
        } else {
            if ( insCorrection == null ) {
                return (snpCorrection == null) ? notFoundResult : new Result(ResultType.SNP_CORRECTED, snpCorrection);
            } else {
                return (snpCorrection != null) ? ambiguousResult : new Result(ResultType.INS_CORRECTED, insCorrection);
            }
        }
    }

    public boolean isInWhitelist(final String seq) {
        return seqs.contains(seq);
    }

    public boolean isInWhitelist(final byte[] bytes, final int start, final int length) {
        return isInWhitelist(new String(Arrays.copyOfRange(bytes, start, start + length)));
    }
}
