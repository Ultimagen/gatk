package org.broadinstitute.hellbender.tools.walkers.pipeline;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.FlowBasedProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.PartialReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

@CommandLineProgramProperties(
        summary = "Split reads file (normally a CRAM) into appropriate read files ready for 10X Cell Ranger invocation",
        oneLineSummary = "Prepare read files for 10X Cell Ranger",
        programGroup = FlowBasedProgramGroup.class
)
@ExperimentalFeature
public class SingleCellPipelineTool extends PartialReadWalker {

    private static final Logger logger = LogManager.getLogger(SingleCellPipelineTool.class);

    // public argument
    @ArgumentCollection
    public SingleCellPipelineToolArgumentCollection args = new SingleCellPipelineToolArgumentCollection();

    // locals
    private SingleCellPipeline pipeline;

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        pipeline.process(read.getName(), read.getBasesNoCopy(), read.getBaseQualitiesNoCopy(), new AttributeProvider() {

            @Override
            public boolean hasAttribute(String attributeName) {
                return read.hasAttribute(attributeName);
            }

            @Override
            public int getAttributeAsInteger(String attributeName) {
                return read.getAttributeAsInteger(attributeName);
            }
        });
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        pipeline = new SingleCellPipeline(args);
        pipeline.onArgsReady();
        pipeline.onStart();
    }

    @Override
    public void closeTool() {
        super.closeTool();
        pipeline.onClose();
    }

    // perform additional argument verification and adjustments, assign defaults
    @Override
    public Object instanceMainPostParseArgs() {

        // validate args
        args.validate();

        return super.instanceMainPostParseArgs();
    }

    @Override
    protected boolean shouldExitEarly(GATKRead read) {

        return pipeline.shouldExitEarly();
    }
}
