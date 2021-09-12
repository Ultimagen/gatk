package org.broadinstitute.hellbender.tools.walkers.annotator.flow;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardMutectAnnotation;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;

@DocumentedFeature(groupName=HelpConstants.DOC_CAT_FLOW_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_FLOW_ANNOTATORS_SUMMARY, summary="Indel Classify Flow Annotation")
public class IndelClassify extends FlowAnnotatorBase implements StandardFlowBasedAnnotation {
    private final static Logger logger = LogManager.getLogger(IndelClassify.class);

    @Override
    public List<String> getKeyNames() {

        return Collections.singletonList(GATKVCFConstants.FLOW_INDEL_CLASSIFY);
    }


}

