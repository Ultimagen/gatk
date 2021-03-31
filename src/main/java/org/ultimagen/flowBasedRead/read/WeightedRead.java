package org.ultimagen.flowBasedRead.read;

public class WeightedRead {
	public final byte[] base; 
	public final byte[] quality; 
	public float weight; 
	WeightedRead( byte[] _base, byte[] _quality, float _weight ){ 
		base=_base;
		quality = _quality;
		weight=_weight; 
	}
	
	public final String toString() {
		String s = new String(base);
		return s + " " + Float.toString(weight); 
	}
}
