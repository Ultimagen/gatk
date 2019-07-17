package org.ultimagenomics.flow_based_read.read;
import org.ultimagenomics.flow_based_read.utils.Variant;

import htsjdk.samtools.SAMRecord;
import java.util.*;
import java.lang.Math;
import java.lang.ArrayIndexOutOfBoundsException; 
// 
public class ReadExpander {
	public List<WeightedRead> expandedReads = null; 
	private List<Variant> variants=null;
	private double nCombinations; 
	public double baseLL = 0;
	public int get_n_haplotypes() {
		return expandedReads.size();
	}
	
	public int get_n_variants() { 
		return variants.size(); 
	}
	
	public final WeightedRead get(int i) {
		if (i < expandedReads.size()) {
			return expandedReads.get(i);
		} else { 
			throw new java.lang.ArrayIndexOutOfBoundsException(String.format("Requested read "
					+ "%d, %d reads available", i, expandedReads.size())); 
		}
	}
	
	
	public ReadExpander(SAMRecord read){
		expandedReads = new ArrayList<WeightedRead>(1);
		if (!read.hasAttribute("LL")) {
			
			expandedReads.add(new WeightedRead(read.getReadBases(), 
					read.getBaseQualities(), 1));
		} else {
			String alt = (String)read.getAttribute("AL"); 
			baseLL = (float)read.getAttribute("LL");

			if ( alt.length() == 0 ){
				expandedReads.add(new WeightedRead(read.getReadBases(), 
						read.getBaseQualities(), 1));
			} else {
				expandedReads.add(new WeightedRead(read.getReadBases(), 
						read.getBaseQualities(), 0));
				
				getVariants(read);
				nCombinations = Math.pow(2, variants.size());
				getAllReadCombinations();
				normalizeWeights();
			}
		}
	}

	
	private void getVariants(SAMRecord read) { 
		String alt = (String) read.getAttribute("AL");
		if( alt.length() > 0 ) {
			String[] mutations = alt.split(";");
			variants = new ArrayList<Variant>(mutations.length);
			for (int j = 0; j<mutations.length; j++) {
				variants.add(j, new Variant(mutations[j]));
			}					
		}
	}
	
	
	private List<Variant> getVariantCombination(int index) {
		if (index >= nCombinations) {
			throw (new ArrayIndexOutOfBoundsException(String.format("Combination %d does not exist", index)));
		}
		List<Variant> result = new ArrayList<Variant>(0);
		for( int cnt = 0; cnt < get_n_variants(); cnt++) {
			if ((index & 1) == 1) { 
				result.add(variants.get(cnt)); 
			}
			index = index >> 1  ; 
		}
		return result; 
	}
	
	private WeightedRead expandRead(List<Variant> variant_combination) {
		
		// length of the new read
		int diff = 0 ;
		float likelihood =  expandedReads.get(0).weight; 
		for (int i = 0 ; i < variant_combination.size(); i++) {
			diff -= (variant_combination.get(i).reference.length() - 
					variant_combination.get(i).alternative.length()); 
			
		}

		byte[] result_seq = new byte[expandedReads.get(0).base.length + diff];
		byte[] result_qual = new byte[expandedReads.get(0).base.length + diff];
		
		int cur_source_pos = 0;
		int cur_dest_pos = 0 ; 
		
		for ( int i = 0; i < variant_combination.size(); i++ ) {

			int update = variant_combination.get(i).position - cur_source_pos;
			if (update >= 0) {
				System.arraycopy(expandedReads.get(0).base,  cur_source_pos, result_seq, 
						cur_dest_pos, update);
				System.arraycopy(expandedReads.get(0).quality,  cur_source_pos, 
						result_qual, cur_dest_pos, update);
			}
			cur_source_pos += update;
			cur_dest_pos += update;
			
			System.arraycopy(variant_combination.get(i).alternative.getBytes(), 0, 
					result_seq, cur_dest_pos, 
					variant_combination.get(i).alternative.length());
			Arrays.fill(result_qual, cur_dest_pos, 
					cur_dest_pos + variant_combination.get(i).alternative.length(), 
					(byte)'I');
			
			cur_source_pos += variant_combination.get(i).reference.length();
			cur_dest_pos += variant_combination.get(i).alternative.length();
			likelihood += variant_combination.get(i).diff;
		}
		
		System.arraycopy(expandedReads.get(0).base,  cur_source_pos, result_seq, 
				cur_dest_pos, expandedReads.get(0).base.length - cur_source_pos);
		System.arraycopy(expandedReads.get(0).quality,  cur_source_pos, result_qual, 
				cur_dest_pos, expandedReads.get(0).base.length-cur_source_pos);
		
		return new WeightedRead(result_seq, result_qual, likelihood );
	}
	
	private void getAllReadCombinations() {
		for (int index = 1 ; index < nCombinations; index++) {
			List<Variant> v = getVariantCombination(index); 
			if (validateVariantCombination(v)) {
				WeightedRead new_read = expandRead(v);
				expandedReads.add(new_read);
			}
		}
		nCombinations = expandedReads.size();
	}
	
	private boolean validateVariantCombination(List<Variant> vc) {
		Set<Integer> posSet = new HashSet<Integer>();
		for (int i = 0 ; i < vc.size(); i++ ) {
			if (posSet.contains(vc.get(i).position))
				return false; 
			posSet.add(vc.get(i).position);
		}
		return true;
	}
	
	private void normalizeWeights() { 
		double total_weight = 0 ; 
		for (int i = 0 ; i < nCombinations; i++ ) { 
			total_weight += Math.pow(2,expandedReads.get(i).weight);
		}
		
		for (int i = 0; i < nCombinations; i++) { 
			expandedReads.get(i).weight = (float)(Math.pow(2, 
					expandedReads.get(i).weight)/total_weight); 
		}
	}
}