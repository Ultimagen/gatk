package org.broadinstitute.hellbender.tools.walkers.annotator.flow;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardMutectAnnotation;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.ultimagen.flowBasedRead.read.FlowBasedRead;

import java.util.*;
import java.util.stream.Collectors;

public abstract class FlowAnnotatorBase extends InfoFieldAnnotation implements StandardMutectAnnotation {
    private final static Logger logger = LogManager.getLogger(FlowAnnotatorBase.class);

    // additional constants
    protected static final String   C_INSERT = "ins";
    protected static final String   C_DELETE = "del";
    protected static final String   C_NA = "NA";
    protected static final String   C_CSS_CS = "cycle-skip";
    protected static final String   C_CSS_PCS = "possible-cycle-skip";
    protected static final String   C_CSS_NS = "non-skip";

    protected static final int      MOTIF_SIZE = 5;
    protected static final int      GC_CONTENT_SIZE = 10;
    protected static final int      BASE_TYPE_COUNT = 4;

    static class LocalAttributes {
        ReferenceContext ref;
        AlleleLikelihoods<GATKRead, Allele> likelihoods;
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
        return annotate(ref, vc, likelihoods, getKeyNames());
    }

    public Map<String, Object> annotate(ReferenceContext ref,
                                        VariantContext vc,
                                        AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                        final List<String> requestedAnnotations) {
        Utils.nonNull(ref);
        Utils.nonNull(vc);

        // some annotators share results
        LocalAttributes         la = new LocalAttributes();
        la.ref = ref;
        la.likelihoods = likelihoods;

        // call annotatotrs
        indelClassify(vc, la);
        if ( requestedAnnotations.contains(GATKVCFConstants.FLOW_HMER_INDEL_LENGTH)
            || requestedAnnotations.contains(GATKVCFConstants.FLOW_HMER_INDEL_NUC)
            || requestedAnnotations.contains(GATKVCFConstants.FLOW_RIGHT_MOTIF) )
            isHmerIndel(vc, la);
        if ( requestedAnnotations.contains(GATKVCFConstants.FLOW_LEFT_MOTIF)
                || requestedAnnotations.contains(GATKVCFConstants.FLOW_RIGHT_MOTIF)
                || requestedAnnotations.contains(GATKVCFConstants.FLOW_CYCLESKIP_STATUS) )
            getMotif(vc, la);
        if ( requestedAnnotations.contains(GATKVCFConstants.FLOW_GC_CONTENT) )
            gcContent(vc, la);
        if ( requestedAnnotations.contains(GATKVCFConstants.FLOW_CYCLESKIP_STATUS) )
            cycleSkip(vc, la);

        // make sure we have a default indel class attibutes
        if ( la.indel == null )
            la.attributes.put(GATKVCFConstants.FLOW_INDEL_CLASSIFY, C_NA);

        // filter map down to requested attributes
        return la.attributes.entrySet().stream()
                .filter(x -> requestedAnnotations.contains(x.getKey()))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    }

    /*
    This function establishes the flow order to be used for manipulating reads in flow space.

    The most natural source for the flow order is are the reads themselves. Alas reads will
    not always be sourced from a bam file with a flow order. In these cases, we can either get it from a
    --flow-order parameter (VariantAnnotator tool) or the default input source (bam)
     */
    private String establishFlowOrder(AlleleLikelihoods<GATKRead, Allele> likelihoods) {

        // extract from a read
        if ( likelihoods != null ) {
            List<GATKRead>  reads = likelihoods.sampleEvidence(0);
            if ( reads.size() > 0 ) {
                GATKRead        read = reads.get(0);
                if ( read instanceof FlowBasedRead ) {
                    return ((FlowBasedRead)read).getFlowOrder();
                } else if ( getFlowOrderForAnnotation() != null )  {
                    establishReadGroupFlowOrder(read.getReadGroup());
                }

            }
        }

        // use global
        return establishReadGroupFlowOrder(null);
    }

    /*
        the flow order might be different for each read group.
        provided flow order can be a list of [group:]flowOrder separated by a comma
        no group: means all/rest
     */
    private String establishReadGroupFlowOrder(String readGroup) {

        for ( String elem : getFlowOrderForAnnotation().split(",") ) {
            String  toks[] = elem.split(":");
            if ( toks.length == 1 )
                return toks[0];
            else if ( toks[0].equals(readGroup) )
                return toks[1];
        }

        throw new RuntimeException("no flow order found for greadGroup " + readGroup);
    }

    // this method should be overriden by annotations requiring a "real" flow order
    // the value returned by this method is a flow order that can be used for computations that work under all flow orders
    protected String getFlowOrderForAnnotation() {
        return "TAGC";
    }

    // "indel_classify" and "indel_length"
    private void indelClassify(final VariantContext vc, final LocalAttributes la) {

        if ( vc.isIndel() ) {

            final List<String>      indelClassify = new LinkedList<>();
            final List<Integer>     indelLength = new LinkedList<>();
            final int               refLength = vc.getReference().length();
            for ( Allele a : vc.getAlleles() ) {
                if ( !a.isReference() ) {
                    indelClassify.add(refLength < a.length() ? C_INSERT : C_DELETE);
                    indelLength.add(Math.abs(refLength - a.length()));
                }
            }
            la.attributes.put(GATKVCFConstants.FLOW_INDEL_CLASSIFY, la.indel = indelClassify);
            la.attributes.put(GATKVCFConstants.FLOW_INDEL_LENGTH, la.indelLength = indelLength);
        }
    }

    /*
    This function determines if the vc is an hmer indel. If so, it marks it as such
     */
    private void isHmerIndel(final VariantContext vc, final LocalAttributes la) {

        // this is (currently) computed only when there is exactly one non reference allele
        if ( vc.isIndel() && la.indel.size() == 1 && vc.getAlleles().size() == 2 ) {

            // establish flow order
            if ( la.flowOrder == null )
                la.flowOrder = establishFlowOrder(la.likelihoods);

            // access alleles
            int         refIndex = vc.getAlleles().get(0).isReference() ? 0 : 1;
            Allele      ref = vc.getAlleles().get(refIndex);
            Allele      alt = vc.getAlleles().get(1 - refIndex);

            // get byte before and after
            byte        before = getReferenceNucleotide(la, vc.getStart() - 1);
            byte[]      after = getReferenceHmerPlus(la, vc.getEnd() + 1, MOTIF_SIZE);

            // build two haplotypes. add byte before and after
            byte[]      refHap = buildHaplotype(before, ref.getBases(), after);
            byte[]      altHap = buildHaplotype(before, alt.getBases(), after);

            // convert to flow space
            int[]       refKey = generateKeyFromSequence(new String(refHap), la.flowOrder, false);
            int[]       altKey = generateKeyFromSequence(new String(altHap), la.flowOrder, false);
            if ( refKey == null || altKey == null )
                throw new GATKException("failed to generate key from reference or alternate sequence");

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

            // check if we've actually encountered a significant different key
            if ( diffIndex < 0 )
                return;
            if ( Math.min(refKey[diffIndex], altKey[diffIndex]) == 0 )
                return;

            // if here, we found the difference.
            byte            nuc = la.flowOrder.getBytes()[diffIndex % la.flowOrder.length()];
            la.hmerIndelLength = refKey[diffIndex];
            la.attributes.put(GATKVCFConstants.FLOW_HMER_INDEL_LENGTH, la.hmerIndelLength);
            la.attributes.put(GATKVCFConstants.FLOW_HMER_INDEL_NUC, Character.toString((char)nuc));

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
        la.attributes.put(GATKVCFConstants.FLOW_LEFT_MOTIF, la.leftMotif = getRefMotif(la, vc.getStart() - MOTIF_SIZE, MOTIF_SIZE));
        if ( vc.isIndel() )
            la.attributes.put(GATKVCFConstants.FLOW_LEFT_MOTIF, la.leftMotif.substring(1) + vc.getReference().getBaseString().substring(0, 1));
        if ( la.rightMotif == null )
            la.rightMotif = getRefMotif(la, vc.getStart() + refLength, MOTIF_SIZE);
        la.attributes.put(GATKVCFConstants.FLOW_RIGHT_MOTIF, la.rightMotif);
    }

    private void gcContent(final VariantContext vc, final LocalAttributes la) {

        int         begin = vc.getStart() - (GC_CONTENT_SIZE / 2);
        String      seq = getRefMotif(la, begin + 1, GC_CONTENT_SIZE);
        int         gcCount = 0;
        for ( byte b : seq.getBytes() ) {
            if ( b == 'G' || b == 'C' ) {
                gcCount++;
            }
        }
        la.attributes.put(GATKVCFConstants.FLOW_GC_CONTENT, (float)gcCount / seq.length());
    }

    private void cycleSkip(final VariantContext vc, final LocalAttributes la) {

        String      css = C_NA;
        if ( !vc.isIndel() && vc.getAlleles().size() == 2 ) {

            // establish flow order
            if ( la.flowOrder == null )
                la.flowOrder = establishFlowOrder(la.likelihoods);

            // access alleles
            int         refIndex = vc.getAlleles().get(0).isReference() ? 0 : 1;
            Allele      ref = vc.getAlleles().get(refIndex);
            Allele      alt = vc.getAlleles().get(1 - refIndex);

            // convert to flow space
            int[]       refKey = generateKeyFromSequence(la.leftMotif + ref.getBaseString() + la.rightMotif, la.flowOrder, true);
            int[]       altKey = generateKeyFromSequence(la.leftMotif + alt.getBaseString() + la.rightMotif, la.flowOrder, true);

            // assign initial css
            css = (refKey.length != altKey.length) ? C_CSS_CS : C_CSS_NS;

            // if same length (NS) then see if it is possible-cycle-skip
            if ( css == C_CSS_NS ) {
                for ( int n = 0 ; n < refKey.length ; n++ ) {
                    if ( (refKey[n] == 0) ^ (altKey[n] == 0) ) {
                        css = C_CSS_PCS;
                        break;
                    }
                }
            }
        }

        la.attributes.put(GATKVCFConstants.FLOW_CYCLESKIP_STATUS, css);
    }

    private int[] generateKeyFromSequence(String sequence, final String flowOrder, boolean ignoreNBases) {

        if ( sequence == null )
            return null;
        if ( ignoreNBases && sequence.indexOf('N') >= 0 )
            sequence = sequence.replace("N", "");

        // allocate maximal key, to be later copied into an array of the exact length.
        int[]       key = new int[sequence.length() * BASE_TYPE_COUNT];
//        
        byte[]      seq = sequence.getBytes();
        byte[]      flow = flowOrder.getBytes();
        int         pos = 0;
        int         keySize = 0;
        for ( int flowPos = 0 ; ; flowPos = (flowPos + 1) % flow.length ) {
            byte    base = flow[flowPos];
            int     hcount = 0;
            for ( int i = pos ; i < seq.length ; i++ ) {
                if ( (seq[i] == 'N') || (seq[i] == '*') )
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
    private byte getReferenceNucleotide(final LocalAttributes la, final int start) {
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
    public void setFlowOrderForTesting(String flowOrder) {

    }
}
