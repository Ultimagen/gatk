package org.broadinstitute.hellbender.tools.walkers.annotator.flow;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AlleleSpecificAnnotationData;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotationData;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.StrandBiasUtils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;

@DocumentedFeature(groupName= HelpConstants.DOC_CAT_FLOW_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_FLOW_ANNOTATORS_SUMMARY, summary="TP based additional error information")
public class AssemblyTPErrorAnnotation extends FlowAnnotatorBase implements StandardFlowBasedAnnotation {

    private final static Logger logger = LogManager.getLogger(AssemblyTPErrorAnnotation.class);
    private static final char ALLELE_SEPARATOR = '|';
    private static final String CALL_SEPARATOR = ",";
    private static DecimalFormat shortDoubleFormat = new DecimalFormat("##0.#####");

    @Override
    public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {

        final LocalContext        localContext = new LocalContext(ref, vc, likelihoods, true);

        if ( localContext.generateAnnotation ) {
            indelClassify(vc, localContext);
            isHmerIndel(vc, localContext);
            variantType(vc, localContext);
            asTpTable(vc, localContext);
        }

        return localContext.asAttributes();
    }

    private void asTpTable(VariantContext vc, LocalContext localContext) {

        if ( localContext.likelihoods != null && localContext.likelihoods.evidenceCount() > 0
                && C_H_MER.equals(localContext.attributes.get(GATKVCFConstants.FLOW_VARIANT_TYPE)) ) {

            // collect reads supporting each allele
            final Map<Allele, List<GATKRead>>     supportingReads = getSupportingReads(vc, localContext);

            // build annotation
            List<String>        annotations = new LinkedList<>();
            for ( Allele allele : supportingReads.keySet() ) {
                final String  annotation = calcAsTpTableAnnotation(vc, localContext, allele, supportingReads.get(allele));
                annotations.add(annotation);
            }
            localContext.attributes.put(GATKVCFConstants.FLOW_AS_TP_TABLE, StringUtils.join(annotations, ALLELE_SEPARATOR));
        } else if ( vc.getAttribute(GATKVCFConstants.FLOW_AS_TP_TABLE) != null ){
            // copy
            localContext.attributes.put(GATKVCFConstants.FLOW_AS_TP_TABLE, vc.getAttribute(GATKVCFConstants.FLOW_AS_TP_TABLE));
        }
    }

    private Map<Allele, List<GATKRead>> getSupportingReads(final VariantContext vc, final LocalContext localContext) {
        final Map<Allele, List<GATKRead>>     supportingReads = new LinkedHashMap<>();
        vc.getAlleles().forEach(allele -> {
            supportingReads.put(allele, new LinkedList<>());
        });
        localContext.likelihoods.samples().forEach(sample -> {
            localContext.likelihoods.bestAllelesBreakingTies(sample).stream()
                    .filter(ba -> ba.isInformative())
                    .forEach(ba -> {
                        supportingReads.get(ba.allele).add(ba.evidence);
                    });
        });
        return supportingReads;
    }

    private boolean isHmerAlleleIsDel(final VariantContext vc, final Allele allele) {
        return allele.length() < vc.getReference().length()
                && vc.getReference().getBaseString().startsWith(allele.getBaseString());
    }

    private String calcAsTpTableAnnotation(final VariantContext vc, final LocalContext localContext,
                                           final Allele allele, final List<GATKRead> reads) {
        final List<String> elems = new LinkedList<>(Arrays.asList(allele.toString(), Integer.toString(reads.size())));
        final Map<String, Integer> counts = new LinkedHashMap<>();
        reads.forEach(read -> {
            final String annotation = calcReadHmerInfo(localContext, vc, allele, read);
            counts.put(annotation, 1 + (counts.containsKey(annotation) ? counts.get(annotation) : 0));
        });
        counts.entrySet().forEach(e -> {
            elems.add(e.getValue() + "x" + e.getKey());
        });
        return StringUtils.join(elems, CALL_SEPARATOR);
    }

    private String calcReadHmerInfo(LocalContext localContext, VariantContext vc, Allele allele, GATKRead read) {
        // convert read into a flow based read
        String              flowOrder = establishFlowOrder(localContext, localContext.likelihoods);
        FlowBasedRead       flowRead = new FlowBasedRead(read, flowOrder, FlowBasedRead.MAX_CLASS, new FlowBasedArgumentCollection());

        // locate flow (index) hmer at the start of the diverging part of the allele
        // assuming the common part of all alleles in the vc is 1 base long
        int                 flow = getReadFlowForLocation(flowRead, vc.getStart() + 1);
        if ( flow < 0 ) {
            return "";
        }
        int                 call = flowRead.getKey()[flow];
        byte                base = flowOrder.getBytes()[flow % flowOrder.length()];

        // get tp for this flow
        final Byte          tp = getReadTPForLocation(flowRead, vc.getStart() + 1);
        if ( tp == null ) {
            return String.format("%d%c", call, base);
        } else {
            return String.format("%d%c/%d", call, base, (byte)tp);
        }
    }

    private int getReadFlowForLocation(FlowBasedRead flowRead, int start) {

        // walk the flow key until offset from start of the read reached
        int     offset = start - flowRead.getStart();
        if ( offset < 0 ) {
            return -1;
        }
        int[]   key = flowRead.getKey();
        for ( int flow = 0 ; flow < key.length ; flow++ ) {
            if ( key[flow] != 0 ) {
                if ( offset < key[flow] ) {
                    return flow;
                }
                offset -= key[flow];
            }
        }

        // if here, not found
        return -1;

    }

    private Byte getReadTPForLocation(FlowBasedRead flowRead, final int start) {

        // walk the flow key until offset from start of the read reached
        int     offset = start - flowRead.getStart();
        if ( offset < 0 ) {
            return null;
        }

        return flowRead.getByteAttributeForBaseIndex(offset, FlowBasedRead.FLOW_MATRIX_TAG_NAME);
    }

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(GATKVCFConstants.FLOW_AS_TP_TABLE);
    }
}
