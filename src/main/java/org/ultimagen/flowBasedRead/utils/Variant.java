package org.ultimagen.flowBasedRead.utils;

public class Variant { 
	public int position;
	public String reference;
	public String alternative;
	public float diff;
	
	public Variant(){ }
	public Variant(String variantTag){
		String[] fields = variantTag.split(",");
		//Note - check if position is zero or one-based
		this.position = Integer.parseInt(fields[0]);
		this.reference = fields[1];
		this.alternative=fields[2]; 
		this.diff = Float.parseFloat(fields[3]);				
	}
	
	public String toString() {
		return String.format("%d, %s -> %s: %f", this.position, 
				this.reference, this.alternative, this.diff);
	}
	
}