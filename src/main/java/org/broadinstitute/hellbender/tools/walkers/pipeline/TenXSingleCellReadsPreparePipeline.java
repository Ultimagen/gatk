package org.broadinstitute.hellbender.tools.walkers.pipeline;

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

import java.util.List;

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
    private String adapter_5p;
    private String adapter_3p;
    private String adapter_middle;

    // temp
    private int readCount;
    private int adapter5pCount;
    private int adapter3pCount;
    private int adapterMiddleCount;


    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {

        // access read
        final byte[]        bases = read.getBasesNoCopy();

        // temp! look for the adapters
        int         adapter5pOfs = AdapterUtils.findAdapter(bases, adapter_5p.getBytes(), args.adapterMinErrorRate, args.adapterMinOverlap);
        int         adapter3pOfs = AdapterUtils.findAdapter(bases, adapter_3p.getBytes(), args.adapterMinErrorRate, args.adapterMinOverlap);
        int         adapterMiddleOfs = AdapterUtils.findAdapter(bases, adapter_middle.getBytes(), args.adapterMinErrorRate, args.adapterMinOverlap);

        // count
        readCount++;
        if ( adapter5pOfs >= 0 ) {
            adapter5pCount++;
        }
        if ( adapter3pOfs >= 0 ) {
            adapter3pCount++;
        }
        if ( adapterMiddleOfs >= 0 ) {
            adapterMiddleCount++;
        }
    }

    @Override
    public void closeTool() {
        super.closeTool();

        logger.info(String.format("counts: %d %d %d %d",
                readCount, adapter5pCount, adapter3pCount, adapterMiddleCount));
    }

    // perform additional argument verification and adjustments, assign defaults
    @Override
    public Object instanceMainPostParseArgs() {

        // validate args
        args.validate();

        // adapters
        if ( (adapter_5p = args.adapter5pOverride) == null ) {
            adapter_5p = (!args.guide || args.libraryDirection == TenXSingleCellArgumentCollection.LibraryDirection.FivePrime)
                    ? "CTACACGACGCTCTTCCGATCT" : "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG";
        }
        if ( (adapter_3p = args.adapter3pOverride) == null ) {
            if ( args.libraryDirection == TenXSingleCellArgumentCollection.LibraryDirection.ThreePrime) {
                adapter_3p = args.guide
                        ? "AGATCGGAAGAGCACACGTCTG" : "CCCATGTACTCTGCGTTGATACCACTGCTT";
            } else {
                adapter_3p = args.guide
                        ? "CTGTCTCTTATACACATCT" : "AGATCGGAAGAGCACACGTCTG";
            }
        }
        if ( (adapter_middle = args.adapterMiddleOverride) == null) {
            if ( args.libraryDirection == TenXSingleCellArgumentCollection.LibraryDirection.FivePrime ) {
                adapter_middle = "^TTTCTTATATGGG";
            } else {
                adapter_middle = args.guide
                        ? "GCTGTTTCCAGCTTAGCTCTTAAAC" : "XTTTTTTTTTTTTTTTTTTTTTTTTT";
            }
        }

        return super.instanceMainPostParseArgs();
    }
}
