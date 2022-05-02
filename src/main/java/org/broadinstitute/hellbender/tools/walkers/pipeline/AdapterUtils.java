package org.broadinstitute.hellbender.tools.walkers.pipeline;

public class AdapterUtils {

    static public int findAdapter(byte[] read, byte[] adapter, double errorRate, int minOverlap) {

        // setup threshold
        int         errorThreshold = (int)(errorRate * adapter.length);

        // trivial implementation to begin with
        for ( int ofs = 0 ; ofs <= read.length - adapter.length ; ofs++ ) {

            // check if an adapter is at this offset
            int     errors = 0;
            for ( int i = 0 ; (i < adapter.length) && (errors < errorThreshold) ; i++ ) {
                if ( !iupacMatch(read[ofs+i], adapter[i]) ) {
                    errors++;
                }
            }

            // found?
            if ( (errors < errorThreshold) && (adapter.length - errors) >= minOverlap ){
                return ofs;
            }
        }

        // if here, did not find
        return -1;
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
