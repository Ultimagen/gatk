package org.broadinstitute.hellbender.tools.walkers.annotator.ultima;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardMutectAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotator;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.ultimagen.flowBasedRead.read.FlowBasedRead;

import java.util.*;
import java.util.stream.Collectors;

@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="UltimaConcordanceAnnotator")
public class UltimaConcordanceAnnotator extends InfoFieldAnnotation implements StandardMutectAnnotation {
    private final static Logger logger = LogManager.getLogger(UltimaConcordanceAnnotator.class);

    // additional constants
    protected static final String   C_INSERT = "ins";
    protected static final String   C_DELETE = "del";
    protected static final String   C_NA = "NA";
    protected static final String   C_CSS_CS = "cycle-skip";
    protected static final String   C_CSS_PCS = "possible-cycle-skip";
    protected static final String   C_CSS_NS = "non-skip";

    protected static final int      MOTIF_SIZE = 5;
    protected static final int      GC_CONTENT_SIZE = 10;

    public static final boolean  addDebugAnnotations = false;

    static class LocalAttributes {
        ReferenceContext ref;
        String      flowOrder;

        List<String> indel;
        List<Integer> indelLength;
        int         hmerIndelLength;
        String      leftMotif;
        String      rightMotif;

        String      refAlleleHmerAndRightMotif;

        Map<String, Object> attributes = new LinkedHashMap<>();
    }

    @Override
    public Map<String, Object> annotate(ReferenceContext ref,
                                        VariantContext vc,
                                        AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        Utils.nonNull(ref);
        Utils.nonNull(vc);

        // some annotators share results
        LocalAttributes         la = new LocalAttributes();
        la.ref = ref;

        // establish flow order
        la.flowOrder = establishFlowOrder(likelihoods);

        // call annotatotrs
        indelClassify(vc, la);
        isHmerIndel(vc, la);
        getMotif(vc, la);
        gcContent(vc, la);
        cycleSkip(vc, la);

        // return attibutes
        if ( la.indel == null )
            la.attributes.put(GATKVCFConstants.ULTIMA_INDEL_CLASSIFY, C_NA);
        if ( addDebugAnnotations ) {
            la.attributes.put(GATKVCFConstants.ULTIMA_DBG_REF, makeDbgRef(ref, vc));
            la.attributes.put(GATKVCFConstants.ULTIMA_DBG_REF_START, ref.getStart());
        }
        return la.attributes;
    }

    private String makeDbgRef(ReferenceContext ref, VariantContext vc) {

        int                 fence = vc.getStart() - ref.getStart();
        StringBuilder       sb = new StringBuilder();

        sb.append(new String(Arrays.copyOfRange(ref.getBases(), 0, fence)));
        sb.append(">");
        sb.append(new String(Arrays.copyOfRange(ref.getBases(), fence, ref.getBases().length)));

        return sb.toString();
    }

    private String establishFlowOrder(AlleleLikelihoods<GATKRead, Allele> likelihoods) {

        // defined?
        if ( VariantAnnotator.flowOrder != null )
            return VariantAnnotator.flowOrder;

            // extract from a read
        if ( likelihoods != null ) {
            List<GATKRead>  reads = likelihoods.sampleEvidence(0);
            if ( reads.size() > 0 ) {
                GATKRead        read = reads.get(0);
                if ( read instanceof FlowBasedRead ) {
                    return ((FlowBasedRead)read).getFlowOrder();
                }
            }
        }

        // if here, it is an error - we have no default for flow order
        throw new RuntimeException("flow-order must be defined. Use --flow-order if running VariantAnnotator");
    }

    @Override
    public List<String> getKeyNames() {

        List<String>        names = new LinkedList<>(Arrays.asList(
                GATKVCFConstants.ULTIMA_INDEL_CLASSIFY, GATKVCFConstants.ULTIMA_INDEL_LENGTH,
                GATKVCFConstants.ULTIMA_HMER_INDEL_LENGTH, GATKVCFConstants.ULTIMA_HMER_INDEL_NUC,
                GATKVCFConstants.ULTIMA_LEFT_MOTIF, GATKVCFConstants.ULTIMA_RIGHT_MOTIF,
                GATKVCFConstants.ULTIMA_GC_CONTENT,
                GATKVCFConstants.ULTIMA_CYCLESKIP_STATUS));

        if ( addDebugAnnotations ) {
            names.addAll(Arrays.asList(
                    GATKVCFConstants.ULTIMA_DBG_REF,
                    GATKVCFConstants.ULTIMA_DBG_REF_START));

        }

        return names;
    }

    // "indel_classify" and "indel_length"
    private void indelClassify(final VariantContext vc, final LocalAttributes la) {

        if ( vc.isIndel() ) {

            /*
            if not x['indel']:
                return None
            elif len(x['ref']) < max([len(y) for y in x['alleles']]):
                return 'ins'
            return 'del'
             */
            final List<String>      indelClassify = new LinkedList<>();
            final List<Integer>     indelLength = new LinkedList<>();
            final int               refLength = vc.getReference().length();
            for ( Allele a : vc.getAlleles() ) {
                if ( !a.isReference() ) {
                    indelClassify.add(refLength < a.length() ? C_INSERT : C_DELETE);
                    indelLength.add(Math.abs(refLength - a.length()));
                }
            }
            la.attributes.put(GATKVCFConstants.ULTIMA_INDEL_CLASSIFY, la.indel = indelClassify);
            la.attributes.put(GATKVCFConstants.ULTIMA_INDEL_LENGTH, la.indelLength = indelLength);
        }
    }

    // "hmer_indel_length" and "hmer_indel_nuc"
    private void isHmerIndel(final VariantContext vc, final LocalAttributes la) {

        // this is (currently) computed only when there is exactly one non reference allele
        if ( vc.isIndel() && la.indel.size() == 1 && vc.getAlleles().size() == 2 ) {

            // access alleles
            int         refIndex = vc.getAlleles().get(0).isReference() ? 0 : 1;
            Allele      ref = vc.getAlleles().get(refIndex);
            Allele      alt = vc.getAlleles().get(1 - refIndex);

            // get byte before and after
            byte        before = getReferenceNucleoid(la, vc.getStart() - 1);
            byte[]      after = getReferenceHmerPlus(la, vc.getEnd() + 1, MOTIF_SIZE);

            // build two haplotypes. add byte before and after
            byte[]      refHap = buildHaplotype(before, ref.getBases(), after);
            byte[]      altHap = buildHaplotype(before, alt.getBases(), after);

            // convert to flow space
            int[]       refKey = generateKeyFromSequence(new String(refHap), la.flowOrder);
            int[]       altKey = generateKeyFromSequence(new String(altHap), la.flowOrder);

            // key must be the same length to begin with
            if ( refKey.length != altKey.length )
                return;

            // key must have only one difference, which should not be between a zero and something
            int     diffIndex = -1;
            int     refBasesCountUpInclHmer = 0;
            for ( int n = 0 ; n < refKey.length ; n++ ) {
                // count ref bases up to and including difference key
                if ( diffIndex < 0 )
                    refBasesCountUpInclHmer += refKey[n];

                // is this the (one) difference key?
                if ( refKey[n] != altKey[n] ) {
                    if ( diffIndex >= 0 )
                        return;
                    else
                        diffIndex = n;
                }
            }
            if ( diffIndex < 0 )
                return;
            if ( Math.min(refKey[diffIndex], altKey[diffIndex]) == 0 )
                return;

            // if here, we found the difference.
            byte            nuc = la.flowOrder.getBytes()[diffIndex % la.flowOrder.length()];
            la.hmerIndelLength = refKey[diffIndex];
            la.attributes.put(GATKVCFConstants.ULTIMA_HMER_INDEL_LENGTH, la.hmerIndelLength);
            la.attributes.put(GATKVCFConstants.ULTIMA_HMER_INDEL_NUC, Character.toString((char)nuc));

            // at this point, we can generate the right motif (for the hmer indel) as we already have the location
            // of the hmer-indel and the bases following it
            la.rightMotif = new String(Arrays.copyOfRange(refHap, refBasesCountUpInclHmer, Math.min(refHap.length, refBasesCountUpInclHmer + MOTIF_SIZE)));
            la.refAlleleHmerAndRightMotif = new String(Arrays.copyOfRange(refHap, 1, Math.min(refHap.length, refBasesCountUpInclHmer + MOTIF_SIZE)));
        }
    }

    private byte[] buildHaplotype(byte before, byte[] bases, byte[] after) {

        byte[]  hap = new byte[1 + bases.length + after.length];

        hap[0] = before;
        System.arraycopy(bases, 0, hap, 1, bases.length);
        System.arraycopy(after, 0, hap, 1 + bases.length, after.length);

        return hap;
    }

    private void getMotif(final VariantContext vc, final LocalAttributes la) {

        // we already did the hard work of building the right motif for hmer-indels. the rest should be simple
        int         refLength = vc.getReference().length();
        la.attributes.put(GATKVCFConstants.ULTIMA_LEFT_MOTIF, la.leftMotif = getRefMotif(la, vc.getStart() - MOTIF_SIZE, MOTIF_SIZE));
        if ( la.rightMotif == null )
            la.attributes.put(GATKVCFConstants.ULTIMA_RIGHT_MOTIF, la.rightMotif = getRefMotif(la, vc.getStart() + refLength, MOTIF_SIZE));
    }

    private void gcContent(final VariantContext vc, final LocalAttributes la) {

        /*
        chrom = faidx[rec['chrom']]
        beg = rec['pos'] - int(size / 2)
        end = beg + size
        seq = chrom[beg:end].seq.upper()
        seqGC = seq.replace('A', '').replace('T', '')
        return float(len(seqGC)) / len(seq)
         */
        int         begin = vc.getStart() - (GC_CONTENT_SIZE / 2);
        String      seq = getRefMotif(la, begin + 1, GC_CONTENT_SIZE);
        int         gcCount = 0;
        for ( byte b : seq.getBytes() ) {
            if ( b == 'G' || b == 'C' ) {
                gcCount++;
            }
        }
        la.attributes.put(GATKVCFConstants.ULTIMA_GC_CONTENT, (float)gcCount / seq.length());
    }

    private void cycleSkip(final VariantContext vc, final LocalAttributes la) {

        /*
        #gt_field is None
        na_pos = df['indel'] | (df['alleles'].apply(len) > 2)
         */
        boolean naIndicator = la.indel != null || vc.getAlleles().size() > 2;

        /*
        df['cycleskip_status'] = "NA"
        snp_pos = ~na_pos
        snps = df.loc[snp_pos].copy()
        left_last = np.array(snps['left_motif']).astype(np.string_)
        right_first = np.array(snps['right_motif']).astype(np.string_)
         */
        String      css = C_NA;
        boolean     snpIndicator = !naIndicator;

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
        if ( la.leftMotif != null && la.rightMotif != null && alt != null ) {
            refSeqs = la.leftMotif + ref + la.rightMotif;
            altSeqs = la.leftMotif + alt + la.rightMotif;
        }

        /*
        ref_encs = [utils.generateKeyFromSequence(
                str(np.char.decode(x)), flow_order) for x in ref_seqs]
        alt_encs = [utils.generateKeyFromSequence(
                str(np.char.decode(x)), flow_order) for x in alt_seqs]
         */
        int[]   refEncs = generateKeyFromSequence(refSeqs, la.flowOrder);
        int[]   altEncs = generateKeyFromSequence(altSeqs, la.flowOrder);

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

            if ( cycleskip )
                css = C_CSS_CS;
            else if ( passCycleskip )
                css = C_CSS_PCS;
            else
                css = C_CSS_NS;
        }

        la.attributes.put(GATKVCFConstants.ULTIMA_CYCLESKIP_STATUS, css);
    }

    private int[] generateKeyFromSequence(final String sequence, final String flowOrder) {

        if ( sequence == null )
            return null;
        int[]         key = new int[sequence.length() * 4];
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
        int         keySize = 0;
        for ( int flowPos = 0 ; ; flowPos = (flowPos + 1) % flow.length ) {
            byte    base = flow[flowPos];
            int     hcount = 0;
            for ( int i = pos ; i < seq.length ; i++ ) {
                if ( seq[i] == 'N')
                    return null;
                else if ( seq[i] == base )
                    hcount++;
                else
                    break;
            }
            if ( pos >= seq.length ) {
                key[keySize++] = hcount;
                break;
            }
            key[keySize++] = hcount;
            pos += hcount;
        }

        return Arrays.copyOfRange(key, 0, keySize);
    }

    // get a single nucleoid from reference
    private byte getReferenceNucleoid(final LocalAttributes la, final int start) {
        int         index = start - la.ref.getWindow().getStart();
        byte[]      bases = la.ref.getBases();
        Utils.validIndex(index, bases.length);
        return bases[index];
    }

    // get an hmer from reference plus a number of additional bases
    private byte[] getReferenceHmerPlus(final LocalAttributes la, final int start, final int additional) {
        int         index = start - la.ref.getWindow().getStart();
        byte[]      bases = la.ref.getBases();
        Utils.validIndex(index, bases.length);

        // get hmer
        StringBuilder sb = new StringBuilder();
        byte          base0 = bases[index++];
        sb.append((char)base0);
        for ( ; index < bases.length && bases[index] == base0 ; index++ )
            sb.append((char)bases[index]);

        // get additional
        for ( int n = 0 ; n < additional && index < bases.length ; n++, index++ )
            sb.append((char)bases[index]);

        return sb.toString().getBytes();
    }
    // get motif from reference
    private String getRefMotif(final LocalAttributes la, final int start, final int length) {
        byte[]      bases = la.ref.getBases();
        int         startIndex = start - la.ref.getWindow().getStart();
        int         endIndex = startIndex + length;
        Utils.validIndex(startIndex, bases.length);
        Utils.validIndex(endIndex-1, bases.length);
        return new String(Arrays.copyOfRange(bases, startIndex, endIndex));
    }

    @VisibleForTesting
    static Map<String, Object> annotateForTesting(final ReferenceContext ref, final VariantContext vc) {

        UltimaConcordanceAnnotator  annotator = new UltimaConcordanceAnnotator();

        return annotator.annotate(ref, vc, null);
    }

}
