package org.broadinstitute.hellbender.tools.walkers.annotator.flow;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.ultimagen.flowBasedRead.read.FlowBasedRead;

import java.util.*;
import java.util.stream.Collectors;

public abstract class FlowAnnotatorBase extends InfoFieldAnnotation {
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
    private static final String     DEFAULT_FLOW_ORDER = "TGCA";

    private List<String>            flowOrder;


    protected class LocalContext {
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

        boolean     notCalculated;

        protected LocalContext(final ReferenceContext ref,
                               final VariantContext vc,
                               final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
            Utils.nonNull(ref);
            Utils.nonNull(vc);

            // some annotators share results
            this.ref = ref;
            this.likelihoods = likelihoods;
        }

        protected Map<String, Object> asAttributes() {

            if ( notCalculated ) {
                return Collections.emptyMap();
            } else {
                return attributes.entrySet().stream()
                        .filter(x -> getKeyNames().contains(x.getKey()))
                        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
            }
        }
    }

    /*
    This function establishes the flow order to be used for manipulating reads in flow space.

    The most natural source for the flow order is are the reads themselves. Alas reads will
    not always be sourced from a bam file with a flow order. In these cases, we can either get it from a
    --flow-order parameter (VariantAnnotator tool) or the default input source (bam)
     */
    private String establishFlowOrder(final LocalContext localContext, final AlleleLikelihoods<GATKRead, Allele> likelihoods) {

        // extract from a read
        if ( likelihoods != null ) {
            final List<GATKRead>  reads = likelihoods.sampleEvidence(0);
            if ( reads.size() > 0 ) {
                GATKRead        read = reads.get(0);
                if ( read instanceof FlowBasedRead ) {
                    return ((FlowBasedRead)read).getFlowOrder();
                } else if ( flowOrder != null )  {
                    establishReadGroupFlowOrder(localContext, read.getReadGroup());
                }
            }
        }

        // use global
        return establishReadGroupFlowOrder(localContext, null);
    }

    /*
        the flow order might be different for each read group.
        provided flow order can be a list of [group:]flowOrder separated by a comma
        no group: means all/rest
     */
    private String establishReadGroupFlowOrder(final LocalContext localContext, final String readGroup) {

        // find flow order for the readGroup
        if ( flowOrder != null ) {
            for (String elem : flowOrder) {
                final String toks[] = elem.split(":");
                if (toks.length == 1) {
                    return toks[0];
                } else if (toks[0].equals(readGroup)) {
                    return toks[1];
                }
            }
        }

        // if here, no flow order was found. may we use a default?
        if ( isActualFlowOrderRequired() ) {
            if ( getNoFlowOrderLogger() != null ) {
                getNoFlowOrderLogger().warn(this.getClass().getSimpleName() + " annotation will not be calculated, no '" + StandardArgumentDefinitions.FLOW_ORDER_FOR_ANNOTATIONS + "' argument provided");
            }
            localContext.notCalculated = true;
        }

        return DEFAULT_FLOW_ORDER;
    }

    protected OneShotLogger getNoFlowOrderLogger() {
        return null;
    }

    protected boolean isActualFlowOrderRequired() {
        return false;
    }

    // "indel_classify" and "indel_length"
    protected void indelClassify(final VariantContext vc, final LocalContext localContext) {

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
            localContext.attributes.put(GATKVCFConstants.FLOW_INDEL_CLASSIFY, localContext.indel = indelClassify);
            localContext.attributes.put(GATKVCFConstants.FLOW_INDEL_LENGTH, localContext.indelLength = indelLength);
        } else {
            localContext.attributes.put(GATKVCFConstants.FLOW_INDEL_CLASSIFY, C_NA);
        }
    }

    /*
    This function determines if the vc is an hmer indel. If so, it marks it as such
     */
    protected void isHmerIndel(final VariantContext vc, final LocalContext localContext) {

        // this is (currently) computed only when there is exactly one non reference allele
        if ( vc.isIndel() && localContext.indel.size() == 1 && vc.getAlleles().size() == 2 ) {

            // establish flow order
            if ( localContext.flowOrder == null ) {
                localContext.flowOrder = establishFlowOrder(localContext, localContext.likelihoods);}


            // access alleles
            final int         refIndex = vc.getAlleles().get(0).isReference() ? 0 : 1;
            final Allele      ref = vc.getAlleles().get(refIndex);
            final Allele      alt = vc.getAlleles().get(1 - refIndex);

            // get byte before and after
            final byte        before = getReferenceNucleotide(localContext, vc.getStart() - 1);
            final byte[]      after = getReferenceHmerPlus(localContext, vc.getEnd() + 1, MOTIF_SIZE);

            // build two haplotypes. add byte before and after
            final byte[]      refHap = buildHaplotype(before, ref.getBases(), after);
            final byte[]      altHap = buildHaplotype(before, alt.getBases(), after);

            // convert to flow space
            final int[]       refKey = generateKeyFromSequence(new String(refHap), localContext.flowOrder, false);
            final int[]       altKey = generateKeyFromSequence(new String(altHap), localContext.flowOrder, false);
            if ( refKey == null || altKey == null ) {
                throw new GATKException("failed to generate key from reference or alternate sequence");
            }

            // key must be the same length to begin with
            if ( refKey.length != altKey.length ) {
                return;
            }

            // key must have only one difference, which should not be between a zero and something
            int     diffIndex = -1;
            int     refBasesCountUpInclHmer = 0;
            for ( int n = 0 ; n < refKey.length ; n++ ) {
                // count ref bases up to and including difference key
                if ( diffIndex < 0 ) {
                    refBasesCountUpInclHmer += refKey[n];
                }

                // is this the (one) difference key?
                if ( refKey[n] != altKey[n] ) {
                    if ( diffIndex >= 0 ) {
                        return;
                    } else {
                        diffIndex = n;
                    }
                }
            }

            // check if we've actually encountered a significant different key
            if ( diffIndex < 0 ) {
                return;
            }
            if ( Math.min(refKey[diffIndex], altKey[diffIndex]) == 0 ) {
                return;
            }

            // if here, we found the difference.
            final byte            nuc = localContext.flowOrder.getBytes()[diffIndex % localContext.flowOrder.length()];
            localContext.hmerIndelLength = refKey[diffIndex];
            localContext.attributes.put(GATKVCFConstants.FLOW_HMER_INDEL_LENGTH, localContext.hmerIndelLength);
            localContext.attributes.put(GATKVCFConstants.FLOW_HMER_INDEL_NUC, Character.toString((char)nuc));

            // at this point, we can generate the right motif (for the hmer indel) as we already have the location
            // of the hmer-indel and the bases following it
            localContext.rightMotif = new String(Arrays.copyOfRange(refHap, refBasesCountUpInclHmer, Math.min(refHap.length, refBasesCountUpInclHmer + MOTIF_SIZE)));
            localContext.refAlleleHmerAndRightMotif = new String(Arrays.copyOfRange(refHap, 1, Math.min(refHap.length, refBasesCountUpInclHmer + MOTIF_SIZE)));
        }
    }

    private byte[] buildHaplotype(final byte before, final byte[] bases, final byte[] after) {

        final byte[]  hap = new byte[1 + bases.length + after.length];

        hap[0] = before;
        System.arraycopy(bases, 0, hap, 1, bases.length);
        System.arraycopy(after, 0, hap, 1 + bases.length, after.length);

        return hap;
    }

    protected void getLeftMotif(final VariantContext vc, final LocalContext localContext) {

        final int         refLength = vc.getReference().length();

        localContext.attributes.put(GATKVCFConstants.FLOW_LEFT_MOTIF, localContext.leftMotif = getRefMotif(localContext, vc.getStart() - MOTIF_SIZE, MOTIF_SIZE));
        if ( vc.isIndel() ) {
            localContext.attributes.put(GATKVCFConstants.FLOW_LEFT_MOTIF, localContext.leftMotif.substring(1) + vc.getReference().getBaseString().substring(0, 1));}
    }

    protected void getRightMotif(final VariantContext vc, final LocalContext localContext) {

        final int         refLength = vc.getReference().length();

        if ( localContext.rightMotif == null ) {
            localContext.rightMotif = getRefMotif(localContext, vc.getStart() + refLength, MOTIF_SIZE);
        }
        localContext.attributes.put(GATKVCFConstants.FLOW_RIGHT_MOTIF, localContext.rightMotif);
    }

    protected void gcContent(final VariantContext vc, final LocalContext localContext) {

        final int         begin = vc.getStart() - (GC_CONTENT_SIZE / 2);
        final String      seq = getRefMotif(localContext, begin + 1, GC_CONTENT_SIZE);
        int         gcCount = 0;
        for ( byte b : seq.getBytes() ) {
            if ( b == 'G' || b == 'C' ) {
                gcCount++;
            }
        }
        localContext.attributes.put(GATKVCFConstants.FLOW_GC_CONTENT, (float)gcCount / seq.length());
    }

    protected void cycleSkip(final VariantContext vc, final LocalContext localContext) {

        // establish flow order
        String      css = C_NA;
        if ( !vc.isIndel() && vc.getAlleles().size() == 2 ) {

            // establish flow order
            if ( localContext.flowOrder == null ) {
                localContext.flowOrder = establishFlowOrder(localContext, localContext.likelihoods);
            }

            // access alleles
            final int         refIndex = vc.getAlleles().get(0).isReference() ? 0 : 1;
            final Allele      ref = vc.getAlleles().get(refIndex);
            final Allele      alt = vc.getAlleles().get(1 - refIndex);

            // convert to flow space
            final int[]       refKey = generateKeyFromSequence(localContext.leftMotif + ref.getBaseString() + localContext.rightMotif, localContext.flowOrder, true);
            final int[]       altKey = generateKeyFromSequence(localContext.leftMotif + alt.getBaseString() + localContext.rightMotif, localContext.flowOrder, true);

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

        localContext.attributes.put(GATKVCFConstants.FLOW_CYCLESKIP_STATUS, css);
    }

    private int[] generateKeyFromSequence(final String inputSequence, final String flowOrder, boolean ignoreNBases) {

        if ( inputSequence == null ) {
            return null;
        }
        final String sequence = (ignoreNBases && inputSequence.indexOf('N') >= 0)
                ? inputSequence.replace("N", "")
                : inputSequence;

        // allocate maximal key, to be later copied into an array of the exact length.
        final int[]       key = new int[sequence.length() * BASE_TYPE_COUNT];
        final byte[]      seq = sequence.getBytes();
        final byte[]      flow = flowOrder.getBytes();
        int         pos = 0;
        int         keySize = 0;
        for ( int flowPos = 0 ; ; flowPos = (flowPos + 1) % flow.length ) {
            byte    base = flow[flowPos];
            int     hcount = 0;
            for ( int i = pos ; i < seq.length ; i++ ) {
                if ( (seq[i] == 'N') || (seq[i] == '*') ) {
                    return null;
                } else if ( seq[i] == base ) {
                    hcount++;
                } else {
                    break;
                }
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
    private byte getReferenceNucleotide(final LocalContext localContext, final int start) {
        final int         index = start - localContext.ref.getWindow().getStart();
        final byte[]      bases = localContext.ref.getBases();
        Utils.validIndex(index, bases.length);
        return bases[index];
    }

    // get an hmer from reference plus a number of additional bases
    private byte[] getReferenceHmerPlus(final LocalContext localContext, final int start, final int additional) {
        int               index = start - localContext.ref.getWindow().getStart();
        final byte[]      bases = localContext.ref.getBases();
        Utils.validIndex(index, bases.length);

        // get hmer
        final StringBuilder sb = new StringBuilder();
        final byte          base0 = bases[index++];
        sb.append((char)base0);
        for ( ; index < bases.length && bases[index] == base0 ; index++ ) {
            sb.append((char) bases[index]);
        }

        // get additional
        for ( int n = 0 ; n < additional && index < bases.length ; n++, index++ ) {
            sb.append((char) bases[index]);
        }

        return sb.toString().getBytes();
    }
    // get motif from reference
    private String getRefMotif(final LocalContext localContext, final int start, final int length) {
        final byte[]      bases = localContext.ref.getBases();
        final int         startIndex = start - localContext.ref.getWindow().getStart();
        final int         endIndex = startIndex + length;
        Utils.validIndex(startIndex, bases.length);
        Utils.validIndex(endIndex-1, bases.length);
        return new String(Arrays.copyOfRange(bases, startIndex, endIndex));
    }

    public void setFlowOrder(final List<String> flowOrder) {
        this.flowOrder = flowOrder;
    }
}
