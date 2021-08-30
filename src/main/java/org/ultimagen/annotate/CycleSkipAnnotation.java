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

        boolean naIndicator = vc.getAttribute(A_INDEL_CLASSIFY) != null || vc.getAlleles().size() > 2;

        String      css = C_CSS_NA;
        boolean     snpIndicator = !naIndicator;
        String      leftLast = vc.getAttributeAsString(A_LEFT_MOTIF, null);
        String      rightFirst = vc.getAttributeAsString(A_RIGHT_MOTIF, null);

        String      ref = vc.getReference().getBaseString();
        String      alt = (snpIndicator && (vc.getAlleles().size() > 1)) ? vc.getAlleles().get(1).getBaseString() : null;

        String      refSeqs = null;
        String      altSeqs = null;
        if ( leftLast != null && rightFirst != null && alt != null ) {
            refSeqs = leftLast + ref + rightFirst;
            altSeqs = leftLast + alt + rightFirst;
        }

        int[]   refEncs = generateKeyFromSequence(refSeqs, flowOrder);
        int[]   altEncs = generateKeyFromSequence(altSeqs, flowOrder);

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
