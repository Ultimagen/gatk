package org.broadinstitute.hellbender.tools.walkers.pipeline;

import java.util.Arrays;

public class AdapterUtils {

    public static final int ADAPTER_NOT_FOUND = -1;

    static public class Adapter {
        final private byte[] adapter;
        final private boolean mustBeAtStart;
        final private boolean mustBeAtEnd;
        final private int errorThreshold;
        final private int minOverlap;
        final private String description;

        public Adapter(final byte[] bases, double errorRate, int minOverlap) {
            this.mustBeAtStart = bases.length > 0 && bases[0] == '^';
            this.mustBeAtEnd = bases.length > 0 && bases[bases.length - 1] == '$';
            this.adapter = Arrays.copyOfRange(bases, mustBeAtStart ? 1 : 0, bases.length - (mustBeAtEnd ? 1 : 0));
            this.errorThreshold = (int)(errorRate * adapter.length);
            this.minOverlap = minOverlap;
            this.description = String.format("%s;max_error_rate=%f;min_overlap=%d", new String(bases), errorRate, minOverlap);
        }

        public int length() {
            return adapter.length;
        }

        public String getDescription() {
            return this.description;
        }

    }

    static public int findAdapter(final byte[] read, final Adapter adapter) {

        // adapter must have some length
        final int adapterLength = adapter.length();
        if ( adapterLength == 0 ) {
            return ADAPTER_NOT_FOUND;
        }


        // trivial implementation to begin with
        int foundOfs = ADAPTER_NOT_FOUND;
        double foundErrorRate = 1.0;
        final int readScanStart = !adapter.mustBeAtEnd ? 0 : read.length - adapterLength;
        final int readScanEnd = read.length - adapterLength;

        for ( int ofs = readScanStart ; ofs <= readScanEnd ; ofs++ ) {

            // check if an adapter is at this offset
            int errors = 0;
            for ( int i = 0 ; (i < adapterLength) && (errors < adapter.errorThreshold) ; i++ ) {
                if ( !iupacMatch(read[ofs+i], adapter.adapter[i]) ) {
                    errors++;
                }
            }

            // exact match found? then return it
            if ( errors == 0 ) {
                return ofs;
            }

            // found?
            if ( (errors < adapter.errorThreshold) && (adapterLength - errors) >= adapter.minOverlap ){
                // update best
                double rate = errors / adapterLength;
                if ( rate < foundErrorRate ) {
                    foundOfs = ofs;
                    foundErrorRate = rate;
                }
            }

            // scanning only at start? if here then it must be the first iteration, let's stop
            if ( adapter.mustBeAtStart ) {
                break;
            }
        }

        // if here, return what we found (or not)
        return foundOfs;
    }

    public static boolean iupacMatch(byte b, byte iupac) {
        if ( iupac == 'A' || iupac == 'C' || iupac == 'G' || iupac == 'T' || iupac == 'U' ) {
            return b == iupac;
        } else if ( iupac == 'X' || iupac == 'N') {
            return b == 'G' || b == 'A' || b == 'T' || b == 'C';
        } else if ( iupac == 'M') {
            return b == 'A' || b == 'C';
        } else if ( iupac == 'R') {
            return b == 'A' || b == 'G';
        } else if ( iupac == 'W') {
            return b == 'A' || b == 'T';
        } else if ( iupac == 'S') {
            return b == 'C' || b == 'G';
        } else if ( iupac == 'Y') {
            return b == 'C' || b == 'T';
        } else if ( iupac == 'K') {
            return b == 'G' || b == 'T';
        } else if ( iupac == 'V') {
            return b == 'A' || b == 'C' || b == 'G';
        } else if ( iupac == 'H') {
            return b == 'A' || b == 'C' || b == 'T';
        } else if ( iupac == 'D') {
            return b == 'A' || b == 'G' || b == 'T';
        } else if ( iupac == 'B') {
            return b == 'C' || b == 'G' || b == 'T';
        } else {
            return false;
        }
    }
}
