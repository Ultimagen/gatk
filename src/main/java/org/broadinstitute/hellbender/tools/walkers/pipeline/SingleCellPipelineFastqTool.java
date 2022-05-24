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
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.*;
import java.util.List;
import java.util.Optional;
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
    protected ProgressMeter progressMeter;


    @Override
    protected Object doWork() {

        // loop on input files
        for ( GATKPath path : input ) {
            try (final BufferedReader reader = new BufferedReader(new InputStreamReader(
                                                new BufferedInputStream(getInputStream(path)))) ) {
                try (final FastqReader fastqReader = new FastqReader(reader)) {

                    for (FastqRecord rec : fastqReader) {
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

                        progressMeter.update(null);
                    }
                }
            } catch (IOException e) {
                throw new GATKException("exception during doWork", e);
            }
        }

        return null;
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
