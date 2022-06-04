package org.broadinstitute.hellbender.tools.walkers.pipeline;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.FlowBasedProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.*;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.zip.GZIPInputStream;

@CommandLineProgramProperties(
        summary = "Split FASTA reads file  into appropriate read files ready for 10X Cell Ranger invocation",
        oneLineSummary = "Prepare read files for 10X Cell Ranger",
        programGroup = FlowBasedProgramGroup.class
)
@ExperimentalFeature
public class SingleCellPipelineFastqTool extends CommandLineProgram {

    private static final Logger logger = LogManager.getLogger(SingleCellPipelineFastqTool.class);

    // public argument
    @ArgumentCollection
    public SingleCellPipelineToolArgumentCollection args = new SingleCellPipelineToolArgumentCollection();

    // fastq input
    @Argument(shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, doc = "input fastq file")
    final private List<GATKPath> input = null;

    // locals
    private SingleCellPipeline pipeline;
    private ProgressMeter progressMeter;
    private BlockingQueue<FastqRecord> fastqRecordQueue;
    private Thread fastqRecordThread;

    @Override
    protected Object doWork() {

        // loop on input files
        for ( GATKPath path : input ) {
            try (final BufferedReader reader = new BufferedReader(new InputStreamReader(
                                                new BufferedInputStream(getInputStream(path)))) ) {
                try (final FastqReader fastqReader = new FastqReader(reader)) {

                    for (FastqRecord rec : fastqReader) {
                        if ( args.multiproc ) {
                            fastqRecordQueue.offer(rec);
                        } else {
                            process(rec);
                        }
                    }
                }
            } catch (IOException e) {
                throw new GATKException("exception during doWork", e);
            }
        }

        // offer end offile and wait for end of processing
        if ( args.multiproc ) {
            fastqRecordQueue.offer(new FastqRecord(null, "", null, ""));
            try {
                fastqRecordThread.join();
            } catch (InterruptedException e) {
                logger.warn("", e);
            }
        }

        return null;
    }

    private void process(final FastqRecord rec) {
        pipeline.lastReason = null;
        pipeline.process(rec.getReadName(), false,  rec.getReadBases(), rec.getBaseQualities(), new AttributeProvider() {
            @Override
            public boolean hasAttribute(String attributeName) {
                return false;
            }

            @Override
            public int getAttributeAsInteger(String attributeName) {
                return 0;
            }
        });
        if ( pipeline.lastReason != null) {
            logger.debug(() -> rec.getReadName() +  "," + pipeline.lastReason);
        }

        progressMeter.update(null);
    }

    private InputStream getInputStream(GATKPath path) throws IOException {
        if ( path.getExtension().isPresent() && path.getExtension().get().equals(".gz") ) {
            return new GZIPInputStream(path.getInputStream());
        } else {
            return path.getInputStream();
        }
    }

    @Override
    protected void onStartup() {
        super.onStartup();

        // progress
        progressMeter = new ProgressMeter();
        progressMeter.start();;

        // create pipeline
        pipeline = new SingleCellPipeline(args);
        pipeline.onArgsReady();
        pipeline.onStart();

        // create queue
        if ( args.multiproc ) {
            fastqRecordQueue = new LinkedBlockingQueue<>(args.queueCapacity);
            (fastqRecordThread = new Thread(() -> {
                try {
                    while (true) {
                        final FastqRecord rec = fastqRecordQueue.take();
                        if (rec.getReadName() == null) {
                            break;
                        } else {
                            process(rec);
                        }
                    }
                } catch (InterruptedException e) {
                    logger.warn("", e);
                }
            })).start();
        }
    }

    @Override
    protected void onShutdown() {
        pipeline.onClose();
        progressMeter.stop();
        super.onStartup();
    }

    // perform additional argument verification and adjustments, assign defaults
    @Override
    public Object instanceMainPostParseArgs() {

        // validate args
        args.validate();

        return super.instanceMainPostParseArgs();
    }
}
