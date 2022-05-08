package org.broadinstitute.hellbender.tools.walkers.pipeline;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.FlowBasedProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.BaseUtils;
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
    private AdapterUtils.AdapterPattern adapter5pPattern;
    private AdapterUtils.AdapterPattern adapter3pPattern;
    private AdapterUtils.AdapterPattern adapterMiddlePattern;

    // fastq output
    private FastqWriterFactory fastqWriterFactory = new FastqWriterFactory();
    private FastqWriter read1Writer;
    private FastqWriter read2Writer;

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {

        if ( args.bypassMode ) {
            read1Writer.write(new FastqRecord(read.getName(), read.getBasesNoCopy(), null, read.getBaseQualitiesNoCopy()));
        } else {
            // access read
            final byte[] bases = read.getBasesNoCopy();

            // temp! look for the adapters
            AdapterUtils.FoundAdapter adapter5p = null;
            AdapterUtils.FoundAdapter adapterMiddle = null;
            AdapterUtils.FoundAdapter adapter3p = null;
            if ( args.bypassMode2 ) {
                final int     length = read.getLength();
                final int     half = length / 2;
                adapter5p = new AdapterUtils.FoundAdapter(0, 1);
                adapterMiddle = new AdapterUtils.FoundAdapter(half, 1);
                adapter3p = new AdapterUtils.FoundAdapter(length - 2, 1);
            } else {
                adapterMiddle = AdapterUtils.findAdapter(bases, adapterMiddlePattern, 0, bases.length, args.returnFirstFoundAdapter);
                if ( adapterMiddle != null ) {
                    adapter5p = AdapterUtils.findAdapter(bases, adapter5pPattern, 0, adapterMiddle.start, args.returnFirstFoundAdapter);
                    if ( adapter5p != null ) {
                        adapter3p = AdapterUtils.findAdapter(bases, adapter3pPattern, adapterMiddle.start + adapterMiddle.length, bases.length, args.returnFirstFoundAdapter);
                    }
                }
            }

            // temp - build 2 reads and write them out
            if (adapter5p != null && adapterMiddle != null && adapter3p != null) {
                if (adapter5p.start >= 0
                        && adapterMiddle.start > (adapter5p.start + adapter5p.length)
                        && adapter3p.start > (adapterMiddle.start + adapterMiddle.length)) {
                    if ( !args.noOutput ) {
                        read1Writer.write(makeFastQRecord(read, adapter5p.start + adapter5p.length, adapterMiddle.start, false));
                        read2Writer.write(makeFastQRecord(read, adapterMiddle.start + adapterMiddle.length, adapter3p.start, args.reverseComplementRead2));
                    }
                }
            }
        }
    }

    private FastqRecord makeFastQRecord(GATKRead read, int startOfs, int endOfs, boolean rc) {

        final int length = endOfs - startOfs;
        byte[] bases = new byte[length];
        final byte[] quals = new byte[length];
        read.copyBases(startOfs, bases, 0, length);
        read.copyBaseQualities(startOfs, quals, 0, length);
        if ( rc ) {
            ArrayUtils.reverse(quals);
            bases = BaseUtils.simpleReverseComplement(bases);
        }

        return new FastqRecord(read.getName(), bases, null, quals);
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        // log adapters
        logger.info("adapter5p: " + adapter5pPattern.getDescription());
        logger.info("adapter3p: " + adapter3pPattern.getDescription());
        logger.info("adapterMiddle: " + adapterMiddlePattern.getDescription());

        // open writers
        File f;
        fastqWriterFactory.setCreateMd5(false);
        fastqWriterFactory.setUseAsyncIo(args.fastqAsyncIO);
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
        String adapter;
        if ( args.no5p3pAdapters ) {
            adapter = "^";
        } else if ( args.adapter5pOverride != null ) {
            adapter = args.adapter5pOverride;
        } else {
            adapter = (!args.guide || args.libraryDirection == TenXSingleCellArgumentCollection.LibraryDirection.FivePrime)
                    ? "CTACACGACGCTCTTCCGATCT" : "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG";
        }
        adapter5pPattern = new AdapterUtils.AdapterPattern(adapter, args.adapterMinErrorRate, args.adapterMinOverlap);

        adapter = null;
        if ( args.no5p3pAdapters ) {
            adapter = "$";
        } else if ( args.adapter3pOverride != null ) {
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
        adapter3pPattern = new AdapterUtils.AdapterPattern(adapter, args.adapterMinErrorRate, args.adapterMinOverlap);

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
        adapterMiddlePattern = new AdapterUtils.AdapterPattern(adapter, args.adapterMinErrorRate, args.adapterMinOverlap);


        return super.instanceMainPostParseArgs();
    }
}
