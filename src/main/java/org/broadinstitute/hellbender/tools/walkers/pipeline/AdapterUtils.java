package org.broadinstitute.hellbender.tools.walkers.pipeline;

import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants;

import java.util.Arrays;

public class AdapterUtils {

    public static final int ADAPTER_NOT_FOUND = -1;

    public static final SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.FASTEST_AVAILABLE);

    static public class AdapterPattern {
        final private byte[] pattern;
        final private boolean mustBeAtStart;
        final private boolean mustBeAtEnd;
        final private int errorThreshold;
        final private int minOverlap;
        final private boolean returnFirstFound;
        final private boolean scanFromEnd;
        final private String description;
        final private String patternString;

        public AdapterPattern(final String pattern, final double errorRate, final int minOverlap,
                              final boolean returnFirstFound, final boolean scanFromEnd) {
            this.mustBeAtStart = pattern.length() > 0 && pattern.charAt(0) == '^';
            this.mustBeAtEnd = pattern.length() > 0 && pattern.charAt(pattern.length() - 1) == '$';
            this.pattern = Arrays.copyOfRange(pattern.getBytes(), mustBeAtStart ? 1 : 0, pattern.length() - (mustBeAtEnd ? 1 : 0));
            this.errorThreshold = (int) (errorRate * this.pattern.length);
            this.minOverlap = minOverlap;
            this.returnFirstFound = returnFirstFound;
            this.scanFromEnd = scanFromEnd;
            this.description = String.format("%s;max_error_rate=%f;min_overlap=%d", pattern, errorRate, minOverlap);
            this.patternString = pattern;
        }

        public int length() {
            return pattern.length;
        }

        public String getDescription() {
            return this.description;
        }

        public byte[] getPattern() {
            return this.pattern;
        }

        public int getErrorThreshold() {
            return errorThreshold;
        }

        public boolean isScanFromEnd() {
            return scanFromEnd;
        }
    }

    static public class FoundAdapter {
        final int start;
        final int length;
        final AdapterPattern adapter;

        public FoundAdapter(final int start, final AdapterPattern adapter) {
            this.start = start;
            this.length = adapter.length();
            this.adapter = adapter;
        }

        public FoundAdapter(final int start, final int length, final AdapterPattern adapter) {
            this.start = start;
            this.length = length;
            this.adapter = adapter;
        }

        public FoundAdapter(final int start, final int length) {
            this(start, length, null);
        }
    }

    static public FoundAdapter findAdapter(final byte[] read, final AdapterPattern adapter, final int start, final int end) {

        // find using simple SNP algorithm
        if ( adapter.scanFromEnd ) {
            return findAdapterWithSNP(read, adapter, start, end);
        } else {
            // if not found, try using alignment
            return findAdapterUsingAligner(read, adapter, start, end);
        }
    }

    static public FoundAdapter findAdapterWithSNP(final byte[] read, final AdapterPattern adapter, final int start, final int end) {

        // adapter must have some length, unless it is a special case
        final int adapterLength = adapter.length();
        if (adapterLength == 0) {
            if (adapter.mustBeAtStart) {
                return new FoundAdapter(0, 0);
            } else if (adapter.mustBeAtEnd) {
                return new FoundAdapter(read.length, 0);
            } else {
                return null;
            }
        }

        // trivial implementation to begin with
        int foundOfs = ADAPTER_NOT_FOUND;
        double foundErrorRate = 1.0;
        final int readScanStart = Math.max(!adapter.mustBeAtEnd ? 0 : read.length - adapterLength, start);
        final int readScanEnd = Math.min(read.length, end) - adapterLength;
        if (readScanStart > readScanEnd) {
            return null;
        }

        // scan from end?
        final int scanIncr;
        final int scanStart;
        final int scanEnd;
        if (!adapter.scanFromEnd) {
            scanIncr = 1;
            scanStart = readScanStart;
            scanEnd = readScanEnd;
        } else {
            scanIncr = -1;
            scanStart = readScanEnd;
            scanEnd = readScanStart;
        }

        // shortcut - a simple forward scan for an exact match
        if ( !adapter.scanFromEnd ) {
            final String readString = new String(read, readScanStart, readScanEnd - readScanStart + adapterLength);
            final int exactMatchOfs = readString.indexOf(adapter.patternString);
            if (exactMatchOfs >= 0) {
                return new FoundAdapter(exactMatchOfs, adapter);
            }
        }

        // scan
        for (int ofs = scanStart; ; ofs += scanIncr) {

            // check if an adapter is at this offset
            int errors = 0;
            for (int i = 0; (i < adapterLength) && (errors <= adapter.errorThreshold); i++) {
                if ((read[ofs + i] != adapter.pattern[i]) && (adapter.pattern[i] != 'X')) {
                    errors++;
                }
            }

            // exact match found? then return it
            if (errors == 0) {
                return new FoundAdapter(ofs, adapter);
            }

            // found?
            if ((errors <= adapter.errorThreshold) && (adapterLength - errors) >= adapter.minOverlap) {
                // update best
                double rate = errors / adapterLength;
                if (rate < foundErrorRate) {
                    foundOfs = ofs;
                    foundErrorRate = rate;
                    if (adapter.returnFirstFound) {
                        break;
                    }
                }
            }

            // scanning only at start? if here then it must be the first iteration, let's stop
            if (adapter.mustBeAtStart) {
                break;
            }

            //  was this the last iteration
            if (ofs == scanEnd) {
                break;
            }
        }

        // if here, return what we found (or not)
        return (foundOfs != ADAPTER_NOT_FOUND) ? new FoundAdapter(foundOfs, adapter) : null;
    }

    public static boolean iupacMatch(byte b, byte iupac) {
        if (iupac == 'A' || iupac == 'C' || iupac == 'G' || iupac == 'T' || iupac == 'U') {
            return b == iupac;
        } else if (iupac == 'X' || iupac == 'N') {
            return b == 'G' || b == 'A' || b == 'T' || b == 'C';
        } else if (iupac == 'M') {
            return b == 'A' || b == 'C';
        } else if (iupac == 'R') {
            return b == 'A' || b == 'G';
        } else if (iupac == 'W') {
            return b == 'A' || b == 'T';
        } else if (iupac == 'S') {
            return b == 'C' || b == 'G';
        } else if (iupac == 'Y') {
            return b == 'C' || b == 'T';
        } else if (iupac == 'K') {
            return b == 'G' || b == 'T';
        } else if (iupac == 'V') {
            return b == 'A' || b == 'C' || b == 'G';
        } else if (iupac == 'H') {
            return b == 'A' || b == 'C' || b == 'T';
        } else if (iupac == 'D') {
            return b == 'A' || b == 'G' || b == 'T';
        } else if (iupac == 'B') {
            return b == 'C' || b == 'G' || b == 'T';
        } else {
            return false;
        }
    }

    static public FoundAdapter findAdapterUsingAligner(final byte[] read, final AdapterPattern adapter, final int start, final int end) {

        final byte[] ref = Arrays.copyOfRange(read, start, end);
        final SmithWatermanAlignment alignment = aligner.align(ref, adapter.getPattern(),
                SmithWatermanAlignmentConstants.STANDARD_NGS, SWOverhangStrategy.SOFTCLIP);

        if ( alignment != null ) {
            // estimate errors (this is not counting SNPs! - also not taking into account min overlap correctly)
            final int   refLength = alignment.getCigar().getReferenceLength();
            if ( refLength >= adapter.minOverlap ) {
                int estimatedErrors = Math.abs(adapter.getPattern().length - refLength);
                if (estimatedErrors < adapter.errorThreshold) {
                    return new FoundAdapter(alignment.getAlignmentOffset(), refLength, adapter);
                }
            }
        }

        // if here, fail
        return null;
    }
}
