package org.broadinstitute.hellbender.tools.walkers.pipeline;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.FlowBasedProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;

@CommandLineProgramProperties(
        summary = "Split reads file (normally a CRAM) into appropriate read files ready for 10X Cell Ranger invocation",
        oneLineSummary = "Prepare read files for 10X Cell Ranger",
        programGroup = FlowBasedProgramGroup.class
)
@ExperimentalFeature
public class TenXSingleCellReadsPreparePipeline extends ReadWalker {

    private static final Logger logger = LogManager.getLogger(TenXSingleCellReadsPreparePipeline.class);

    // public argument
    @ArgumentCollection
    public TenXSingleCellArgumentCollection args = new TenXSingleCellArgumentCollection();

    // locals
    private AdapterUtils.Adapter adapter5p;
    private AdapterUtils.Adapter adapter3p;
    private AdapterUtils.Adapter adapterMiddle;

    // fastq output
    private FastqWriterFactory fastqWriterFactory = new FastqWriterFactory();
    private FastqWriter read1Writer;
    private FastqWriter read2Writer;

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {

        // access read
        final byte[]        bases = read.getBasesNoCopy();

        // temp! look for the adapters
        int         adapter5pOfs = AdapterUtils.findAdapter(bases, adapter5p);
        int         adapter3pOfs = AdapterUtils.findAdapter(bases, adapter3p);
        int         adapterMiddleOfs = AdapterUtils.findAdapter(bases, adapterMiddle);

        // temp - build 2 reads and write them out
        if ( adapter5pOfs >= 0
                && adapterMiddleOfs > (adapter5pOfs + adapter5p.length())
                && adapter3pOfs > (adapterMiddleOfs + adapterMiddle.length()) ) {
            read1Writer.write(makeFastQRecord(read, adapter5pOfs + adapter5p.length(), adapterMiddleOfs));
            read2Writer.write(makeFastQRecord(read, adapterMiddleOfs + adapterMiddle.length(), adapter3pOfs));
        }
    }

    private FastqRecord makeFastQRecord(GATKRead read, int startOfs, int endOfs) {

        final int length = endOfs - startOfs;
        final byte[] bases = new byte[length];
        final byte[] quals = new byte[length];
        read.copyBases(startOfs, bases, 0, length);
        read.copyBaseQualities(startOfs, quals, 0, length);

        return new FastqRecord(read.getName(), bases, null, quals);
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        // log adapters
        logger.info("adapter5p: " + adapter5p.getDescription());
        logger.info("adapter3p: " + adapter3p.getDescription());
        logger.info("adapterMiddle: " + adapterMiddle.getDescription());

        // open writers
        File f;
        read1Writer = fastqWriterFactory.newWriter(f = buildFastQReedsOutputFile(1));
        logger.info("read1 output: " + f);
        read2Writer = fastqWriterFactory.newWriter(f = buildFastQReedsOutputFile(2));
        logger.info("read2 output: " + f);

    }

    private File buildFastQReedsOutputFile(int i) {
        return new File(args.baseFilename + "_read" + i + ".fastq");
    }

    @Override
    public void closeTool() {
        super.closeTool();

        // close writers
        if ( read1Writer != null ) {
            read1Writer.close();;
        }
        if ( read2Writer != null ) {
            read2Writer.close();;
        }
    }

    // perform additional argument verification and adjustments, assign defaults
    @Override
    public Object instanceMainPostParseArgs() {

        // validate args
        args.validate();

        // adapters
        String adapter = null;
        if ( args.adapter5pOverride != null ) {
            adapter = args.adapter5pOverride;
        } else {
            adapter = (!args.guide || args.libraryDirection == TenXSingleCellArgumentCollection.LibraryDirection.FivePrime)
                    ? "CTACACGACGCTCTTCCGATCT" : "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG";
        }
        adapter5p = new AdapterUtils.Adapter(adapter.getBytes(), args.adapterMinErrorRate, args.adapterMinOverlap);

        adapter = null;
        if ( args.adapter3pOverride != null ) {
            adapter = args.adapter3pOverride;
        } else {
            if ( args.libraryDirection == TenXSingleCellArgumentCollection.LibraryDirection.ThreePrime) {
                adapter = args.guide
                        ? "AGATCGGAAGAGCACACGTCTG" : "CCCATGTACTCTGCGTTGATACCACTGCTT";
            } else {
                adapter = args.guide
                        ? "CTGTCTCTTATACACATCT" : "AGATCGGAAGAGCACACGTCTG";
            }
        }
        adapter3p = new AdapterUtils.Adapter(adapter.getBytes(), args.adapterMinErrorRate, args.adapterMinOverlap);

        adapter = null;
        if ( args.adapterMiddleOverride != null ) {
            adapter = args.adapterMiddleOverride;
        } else {
            if ( args.libraryDirection == TenXSingleCellArgumentCollection.LibraryDirection.FivePrime ) {
                adapter = "^TTTCTTATATGGG";
            } else {
                adapter = args.guide
                        ? "GCTGTTTCCAGCTTAGCTCTTAAAC" : "XTTTTTTTTTTTTTTTTTTTTTTTTT";
            }
        }
        adapterMiddle = new AdapterUtils.Adapter(adapter.getBytes(), args.adapterMinErrorRate, args.adapterMinOverlap);


        return super.instanceMainPostParseArgs();
    }
}
