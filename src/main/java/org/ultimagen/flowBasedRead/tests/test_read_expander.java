package org.ultimagen.flowBasedRead.tests;

import java.io.*;

import htsjdk.samtools.*;
import org.ultimagen.flowBasedRead.read.ReadExpander;

public class test_read_expander {
	public static void main(String args[]) {
		
		System.out.println("Input file: " + args[0]);
		final SamReader reader = SamReaderFactory.makeDefault().open(new File(args[0]));
		int cnt = 0 ;
		
		for (final SAMRecord samRecord : reader ) {
			//System.out.println("***************************");
			if (samRecord.hasAttribute("LL")){
				ReadExpander rexp = new ReadExpander(samRecord);
				for ( int i = 0 ; i < rexp.get_n_haplotypes() ; i++ ) {
					
					;//System.out.println(rexp.get(i).toString());
				}
			}
			
		}
	}

}
