package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrReferenceAnalyzer;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

/**
 * This is the genotyper model that genotypes homopolymer length.
 *
 *
 */

public class HmerGenotypesModel implements GenotypingModel {


    @Override
    public <A extends Allele> GenotypingLikelihoods<A> calculateLikelihoods(AlleleList<A> genotypingAlleles, GenotypingData<A> data,
                                                                            byte[] paddedReference, int offsetForRefIntoEvent,
                                                                            DragstrReferenceAnalyzer dragstrs) {
        return null;
    }

    static public boolean isEligibleHomopolymerIndel(final VariantContext vc, final int loc,
                                                     final DragstrReferenceAnalyzer dragstrs, final int hpolIndelThreshold) {
        if (vc==null) {
            return false;
        }

        final int period = dragstrs.period(loc);
        final int repeats = dragstrs.repeatLength(loc);
        final byte ru = dragstrs.repeatUnit(loc)[0];
        if ((period == 1) && (repeats >= hpolIndelThreshold)){
            if (!vc.isIndel() || !vc.getAlleles().stream().allMatch(a -> isHmerIndel(a,ru))) {
                return false;
            }
        } else {
            return false;
        }
        return true;
    }


    static protected boolean isHmerIndel(final Allele al, final byte hmer_base){
        for (int i = 1; i< al.length(); i++){
            if (al.getBases()[i] != hmer_base){
                return false;
            }
        }
        return true;
    }

    public static boolean isRefBiasExpected(final Event allele, final Haplotype refHaplotype, final int biasedIndelThreshold) {
        Haplotype altHaplotype = refHaplotype.insertAllele(
                allele.refAllele(),
                allele.altAllele(),
                allele.getStart());
        Pair<Integer, Integer> indel_length_and_type = BaseUtils.equalUpToHmerChangeGetLength(refHaplotype.getBases(), altHaplotype.getBases());
        if (indel_length_and_type.getLeft() == 0) {
            return false;
        } else if (indel_length_and_type.getLeft() == -1) {
            return false;
        } else if (indel_length_and_type.getRight() > biasedIndelThreshold){
            return true;
        }
        return true;
    }

}
