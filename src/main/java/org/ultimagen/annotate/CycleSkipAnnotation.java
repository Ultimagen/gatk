package org.ultimagen.annotate;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CycleSkipAnnotation extends AnnotatorBase {

    private String      flowOrder = "TACG";

    public CycleSkipAnnotation(final ReferenceSequenceFile referenceSequenceFile) {
        super(referenceSequenceFile);
    }

    @Override
    public void annotate(final VariantContext vc) {

        /*
        #gt_field is None
        na_pos = df['indel'] | (df['alleles'].apply(len) > 2)
         */
        boolean naIndicator = vc.getAttribute(A_INDEL_CLASSIFY) != null || vc.getAlleles().size() > 2;

        /*
        df['cycleskip_status'] = "NA"
        snp_pos = ~na_pos
        snps = df.loc[snp_pos].copy()
        left_last = np.array(snps['left_motif']).astype(np.string_)
        right_first = np.array(snps['right_motif']).astype(np.string_)
         */
        String      css = C_CSS_NA;
        boolean     snpIndicator = !naIndicator;
        String      leftLast = vc.getAttributeAsString(A_LEFT_MOTIF, null);
        String      rightFirst = vc.getAttributeAsString(A_RIGHT_MOTIF, null);

        /*
        ref = np.array(snps['ref']).astype(np.string_)
        alt = np.array(snps['alleles'].apply(
            lambda x: x[1] if len(x) > 1 else None)).astype(np.string_)
         */
        String      ref = vc.getReference().getBaseString();
        String      alt = (snpIndicator && (vc.getAlleles().size() > 1)) ? vc.getAlleles().get(1).getBaseString() : null;

        /*
        ref_seqs = np.char.add(np.char.add(left_last, ref), right_first)
        alt_seqs = np.char.add(np.char.add(left_last, alt), right_first)
         */
        String      refSeqs = null;
        String      altSeqs = null;
        if ( leftLast != null && rightFirst != null && alt != null ) {
            refSeqs = leftLast + ref + rightFirst;
            altSeqs = leftLast + alt + rightFirst;
        }

        /*
        ref_encs = [utils.generateKeyFromSequence(
                str(np.char.decode(x)), flow_order) for x in ref_seqs]
        alt_encs = [utils.generateKeyFromSequence(
                str(np.char.decode(x)), flow_order) for x in alt_seqs]
         */
        int[]   refEncs = generateKeyFromSequence(refSeqs, flowOrder);
        int[]   altEncs = generateKeyFromSequence(altSeqs, flowOrder);

        /*
        cycleskip = np.array([x for x in range(len(ref_encs))
                          if len(ref_encs[x]) != len(alt_encs[x])])
        poss_cycleskip = [x for x in range(len(ref_encs)) if len(ref_encs[x]) == len(alt_encs[x])
                          and (np.any(ref_encs[x][ref_encs[x] - alt_encs[x] != 0] == 0) or
                           np.any(alt_encs[x][ref_encs[x] - alt_encs[x] != 0] == 0))]
        s = set(np.concatenate((cycleskip, poss_cycleskip)))
        non_cycleskip = [x for x in range(len(ref_encs)) if x not in s]
         */
        if ( altEncs != null ) {
            boolean cycleskip = refEncs.length != altEncs.length;
            boolean passCycleskip = !cycleskip && !Arrays.equals(refEncs, altEncs);
            boolean s = cycleskip || passCycleskip;
            boolean nonCycleSkip = !s;

            if ( cycleskip )
                css = C_CSS_CS;
            else if ( passCycleskip )
                css = C_CSS_PCS;
            else if ( nonCycleSkip )
                css = C_CSS_NS;
        }

        annotate(vc, A_CYCLESKIP_STATUS, css);
    }

    private int[] generateKeyFromSequence(final String sequence, final String flowOrder) {

        if ( sequence == null )
            return null;
        List<Integer>       key = new ArrayList<>(sequence.length() * (flowOrder.length() - 1));

        /*
        key = []
        pos = 0
        for base in flow:
            hcount = 0
            for i in range(pos, len(sequence)):
                if sequence[i] == base:
                    hcount+=1
                else:
                    break
            else:
                key.append(hcount)
                break # end of sequence
            key.append(hcount)
            pos += hcount
         */
        byte[]      seq = sequence.getBytes();
        byte[]      flow = flowOrder.getBytes();
        int         pos = 0;
        for ( int flowPos = 0 ; ; flowPos = (flowPos + 1) % flow.length ) {
            byte    base = flow[flowPos];
            int     hcount = 0;
            for ( int i = pos ; i < seq.length ; i++ ) {
                if ( seq[i] == base )
                    hcount++;
                else
                    break;
            }
            if ( pos >= seq.length ) {
                key.add(hcount);
                break;
            }
            key.add(hcount);
            pos += hcount;
        }

        return key.stream().mapToInt(i->i).toArray();
    }
}
