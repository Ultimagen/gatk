package org.ultimagen.flowBasedRead.read;

import java.util.ArrayList;

public class FlowBasedKeyCodec {

    /**
     * Converts base space sequence to flow space
     * @param bases base space sequence
     * @param flowOrder flow order
     * @param clipping maximal flow value to output (at most 127).
     * @return Array of flow values
     */
    static protected byte[] base2key(final byte[] bases, final String flowOrder, int clipping) {

        clipping = Math.min(clipping, 127);

        final int[]       intKeys = base2key(bases, flowOrder);
        final byte[]      byteKeys = new byte[intKeys.length];
        int         i = 0;

        for ( final int intKey : intKeys ) {
            byteKeys[i++] = (byte)((intKey < clipping) ? intKey : clipping);
        }

        return byteKeys;
    }

    /**
     * Converts base space sequence to flow space
     * @param bases base space sequence
     * @param flowOrder flow order
     * @return Array of flow values
     */

    static public int[] base2key(final byte[] bases, final String flowOrder){

        final ArrayList<Integer> result = new ArrayList<>();
        final byte[] flowOrderBytes = flowOrder.getBytes();
        int loc = 0;
        int flowNumber = 0 ;
        final int period = flowOrderBytes.length;
        while ( loc < bases.length ) {
            final byte flowBase = flowOrderBytes[flowNumber%period];
            if ((bases[loc]!=flowBase) && ( bases[loc]!= FlowBasedRead.N_ASCII)) {
                result.add(0);
            } else {
                int count = 0;
                while ( ( loc < bases.length) && ((bases[loc]==flowBase) || (bases[loc]== FlowBasedRead.N_ASCII)) ){
                    loc++;
                    count ++;
                }
                result.add(count);

            }
            flowNumber++;
        }
        final int[] ret = new int[result.size()];
        for (int i = 0; i < result.size(); i++) {
            ret[i] = result.get(i);
        }
        return ret;
    }

    /**
     * For every flow of the key output the index of the last base that was output prior to this flow
     * @param key
     * @return array
     */
    static protected int[] getKey2Base(final byte[] key){
        final int[] result = new int[key.length];
        result[0] = -1;
        for (int i = 1 ; i < result.length; i++) {
            result[i] = result[i-1] + key[i-1];
        }
        return result;
    }

    /**
     * For every flow of the key output the index of the last base that was output prior to this flow
     * @param key
     * @return array
     */
    static protected int[] getKey2Base(final int[] key) {
        final int[] result = new int[key.length];
        result[0] = -1;
        for (int i = 1; i < result.length; i++) {
            result[i] = result[i - 1] + key[i - 1];
        }
        return result;

    }

    /**
     * For every flow of the key output the nucleotide that is being read for this flow
     * @param flowOrder
     * @param expectedLength the length of the key (key is not provided)
     * @return array of bases
     */

    static protected byte[] getFlow2Base(final String flowOrder, final int expectedLength) {
        final byte[] result = new byte[expectedLength] ;
        for ( int i = 0; i < result.length; i++ ) {
            result[i] = (byte)flowOrder.charAt(i%flowOrder.length());
        }
        return result;
    }

    /**
     * Prints the key as character-encoded string
     * @param bytes  key array
     * @return encoding
     */
    static public String keyAsString(final byte[] bytes)
    {
        final StringBuilder   sb = new StringBuilder();

        for ( final byte b : bytes )
            sb.append((char)((b < 10) ? ('0' + b) : ('A' + b - 10)));

        return sb.toString();
    }

    /**
     * Prints the key as character-encoded string
     * @param ints (key array)
     * @return encoded string
     */
    static public String keyAsString(final int[] ints)
    {
        final StringBuilder   sb = new StringBuilder();

        for ( final int i : ints )
            sb.append((char)((i < 10) ? ('0' + i) : ('A' + i - 10)));

        return sb.toString();
    }
}
