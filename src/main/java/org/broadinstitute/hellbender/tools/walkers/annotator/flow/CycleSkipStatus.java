package org.broadinstitute.hellbender.tools.walkers.annotator.flow;

import com.google.common.annotations.VisibleForTesting;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardMutectAnnotation;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Collections;
import java.util.List;

@DocumentedFeature(groupName=HelpConstants.DOC_CAT_FLOW_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_FLOW_ANNOTATORS_SUMMARY, summary="Cycle Skip Status Flow Annotation")
public class CycleSkipStatus extends FlowAnnotatorBase implements StandardFlowBasedAnnotation {
    private final static Logger logger = LogManager.getLogger(CycleSkipStatus.class);

    @Argument(fullName = "flow-order-for-cycle-skip",  doc = "flow order used for this annotations. [readGroup:]flowOrder,...", optional = true)
    private String flowOrder;

    @Override
    public List<String> getKeyNames() {

        return Collections.singletonList(GATKVCFConstants.FLOW_CYCLESKIP_STATUS);
    }

    @Override
    protected String getFlowOrderForAnnotation() {
        if ( flowOrder != null )
            return flowOrder;

        // if here, it is an error - we have no default for flow order
        throw new RuntimeException("flow order must be defined. Use --flow-order-for-annotations parameter");
    }

    @VisibleForTesting
    @Override
    public void setFlowOrderForTesting(String flowOrder) {
        this.flowOrder = flowOrder;
    }

}

