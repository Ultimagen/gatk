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
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.Collections;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

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
    private BlockingQueue<GATKRead> readQueue;
    private Thread readThread;

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        if ( args.multiproc ) {
            readQueue.offer(read);
        } else {
            process(read);
        }
    }

    private void process(final GATKRead read) {

        pipeline.process(read.getName(), read.isReverseStrand(), read.getBasesNoCopy(), read.getBaseQualitiesNoCopy(), new AttributeProvider() {

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
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.emptyList();
    }

    @Override
    public void traverse() {
        super.traverse();

        // wait for processing to be over
        readQueue.offer(new SAMRecordToGATKReadAdapter(null));
        try {
            readThread.join();
        } catch (InterruptedException e) {
            logger.warn("", e);
        }
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        pipeline = new SingleCellPipeline(args);
        pipeline.onArgsReady();
        pipeline.onStart();

        // create queue
        if ( args.multiproc ) {
            readQueue = new LinkedBlockingQueue<>(args.queueCapacity);
            (readThread = new Thread(() -> {
                try {
                    while (true) {
                        final GATKRead read = readQueue.take();
                        if ((read instanceof SAMRecordToGATKReadAdapter) && ((SAMRecordToGATKReadAdapter)read).getEncapsulatedSamRecord() == null) {
                            break;
                        } else {
                            process(read);
                        }
                    }
                } catch (InterruptedException e) {
                    logger.warn("", e);
                }
            })).start();
        }
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
