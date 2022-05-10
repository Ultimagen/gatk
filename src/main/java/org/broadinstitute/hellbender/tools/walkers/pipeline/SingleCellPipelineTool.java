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
import org.broadinstitute.hellbender.engine.PartialReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;

@CommandLineProgramProperties(
        summary = "Split reads file (normally a CRAM) into appropriate read files ready for 10X Cell Ranger invocation",
        oneLineSummary = "Prepare read files for 10X Cell Ranger",
        programGroup = FlowBasedProgramGroup.class
)
@ExperimentalFeature
public class SingleCellPipelineTool extends PartialReadWalker {

    private static final Logger logger = LogManager.getLogger(SingleCellPipelineTool.class);
    public static final int CBC_LENGTH = 16;
    public static final int CBC_UMI_MAX_LENGTH = 390;
    public static final int CDNA_MAX_LENGTH = 315;
    public static final String REPORT_JSON_FILEnAME_SUFFIX = "_report.json";

    // public argument
    @ArgumentCollection
    public SingleCellPipelineToolArgumentCollection args = new SingleCellPipelineToolArgumentCollection();

    // adapters
    private AdapterUtils.AdapterPattern adapter5pPattern;
    private AdapterUtils.AdapterPattern adapter3pPattern;
    private AdapterUtils.AdapterPattern adapterMiddlePattern;

    // fastq output
    private FastqWriterFactory fastqWriterFactory = new FastqWriterFactory();
    private FastqWriter read1Writer;
    private FastqWriter read2Writer;

    // other locals
    private SingleCellPipelineToolStatistics stats = new SingleCellPipelineToolStatistics();

    private int umiLengthOrig;
    private int cbcUmiLength;
    private int cbcUmiLengthOrig;

    static class FoundAdapters {
        AdapterUtils.FoundAdapter adapter5p;
        AdapterUtils.FoundAdapter adapterMiddle;
        AdapterUtils.FoundAdapter adapter3p;

        FoundAdapters(final AdapterUtils.FoundAdapter adapter5p, final AdapterUtils.FoundAdapter adapterMiddle, final AdapterUtils.FoundAdapter adapter3p) {
            this.adapter5p = adapter5p;
            this.adapterMiddle = adapterMiddle;
            this.adapter3p = adapter3p;
        }
    }

    @Override
    public void apply(GATKRead read, ReferenceContext referenceContext, FeatureContext featureContext) {

        stats.readsIn++;
        stats.bpIn += read.getLength();

        if ( args.bypassMode ) {
            read1Writer.write(new FastqRecord(read.getName(), read.getBasesNoCopy(), null, read.getBaseQualitiesNoCopy()));
        } else {
            // rsq threshold?
            if ( args.rsqThreshold > 0 && read.hasAttribute("RQ")
                && read.getAttributeAsInteger("RQ") < args.rsqThreshold ) {
                stats.rqDropped++;
                return;
            }

            // access read
            final byte[] bases = read.getBasesNoCopy();
            final int basesTrimmedLength = findTrimmedLength(bases, read.getBaseQualitiesNoCopy(), args.qualityCutoff);
            stats.bpCutoff += (bases.length - basesTrimmedLength);

            // find adapters
            FoundAdapters foundAdapters = findAdapters(read, bases, basesTrimmedLength);

            if ( foundAdapters != null ) {

                // find and limit read1 (cbc_umi)
                final int read1Start = foundAdapters.adapter5p.start + foundAdapters.adapter5p.length;
                final int read1Length = Math.min(foundAdapters.adapterMiddle.start - read1Start, cbcUmiLengthOrig);
                final boolean read1Valid = read1Length == cbcUmiLengthOrig;
                stats.read1TooShortDropped += (read1Valid ? 0 : 1);

                // find and limit read1 (cdna)
                int read2Start = foundAdapters.adapterMiddle.start + foundAdapters.adapterMiddle.length;
                int read2Length = Math.min(foundAdapters.adapter3p.start - read2Start, CDNA_MAX_LENGTH);
                boolean read2Valid = read2Length >= args.minCdnaLength;
                if ( args.cdnaFirstBasesToClip != 0 ) {
                    read2Start += args.cdnaFirstBasesToClip;
                    read2Length -= args.cdnaFirstBasesToClip;
                }
                if ( args.cdnaTrimmingLength != 0 ) {
                    read2Length = Math.min(read2Length, args.cdnaTrimmingLength);
                }
                if ( read2Length <= 0 ) {
                    read2Valid = false;
                }
                stats.read2TooShortDropped += (read2Valid ? 0 : 1);

                // reads valid?
                if (read1Valid && read2Valid && !args.noOutput) {

                    if (args.logAdapters == SingleCellPipelineToolArgumentCollection.LogAdapters.Output) {
                        logAdapters(read, foundAdapters);
                    }

                    read1Writer.write(makeFastQRecord(read, read1Start, read1Start + read1Length, false, true));
                    read2Writer.write(makeFastQRecord(read, read2Start, read2Start + read2Length, args.reverseComplementRead2, false));

                    stats.readsOut++;
                    stats.bpOut += (read1Length + read2Length);
                }
            }
        }
    }

    private FoundAdapters findAdapters(GATKRead read, byte[] bases, int basesTrimmedLength) {

        FoundAdapters result = null;

        // temp! look for the adapters
        final AdapterUtils.FoundAdapter adapter5p = AdapterUtils.findAdapter(bases, adapter5pPattern, 0, basesTrimmedLength);
        final AdapterUtils.FoundAdapter adapter3p = AdapterUtils.findAdapter(bases, adapter3pPattern, 0, basesTrimmedLength);
        final AdapterUtils.FoundAdapter adapterMiddle = AdapterUtils.findAdapter(bases, adapterMiddlePattern, 0, basesTrimmedLength);
        if ( adapter5p != null && adapter3p != null && adapterMiddle != null ) {
            result = new FoundAdapters(adapter5p, adapterMiddle, adapter3p);
        }

        // stats
        stats.adapter5p += (adapter5p != null ? 1 : 0);
        stats.adapterMiddle += (adapterMiddle != null ? 1 : 0);
        stats.adapter3p += (adapter3p != null ? 1 : 0);

        // log adapters?
        if ( args.logAdapters == SingleCellPipelineToolArgumentCollection.LogAdapters.Input ) {
            logAdapters(read,
                    new AdapterUtils.FoundAdapter[] {adapter5p, adapterMiddle, adapter3p});
        }

        return result;

    }

    private int findTrimmedLength(byte[] bases, byte[] quals, int qualityCutoff) {
        int     length = bases.length;
        if ( qualityCutoff != 0 ) {
            while ( length > 0 && quals[length-1] < qualityCutoff ) {
                length--;
            }
        }
        return length;
    }

    private FastqRecord makeFastQRecord(GATKRead read, int startOfs, int endOfs, boolean rc, boolean maskLast) {

        final int length = endOfs - startOfs;
        byte[] bases = new byte[length];
        final byte[] quals = new byte[length];
        read.copyBases(startOfs, bases, 0, length);
        read.copyBaseQualities(startOfs, quals, 0, length);
        if ( maskLast && args.cbcUmiMaskLastBytes > 0 ) {
            for ( int i = 0 ; i < args.cbcUmiMaskLastBytes ; i++ ) {
                bases[bases.length - i - 1] = 'T';
                quals[quals.length - i - 1] = 'I';
            }
        }
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

        // write stats
        stats.writeJson(new File(args.baseFilename + REPORT_JSON_FILEnAME_SUFFIX));
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
            adapter = (!args.guide || args.libraryDirection == SingleCellPipelineToolArgumentCollection.LibraryDirection.FivePrime)
                    ? "CTACACGACGCTCTTCCGATCT" : "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG";
        }
        adapter5pPattern = new AdapterUtils.AdapterPattern(adapter, args.adapterMinErrorRate, args.adapterMinOverlap, args.returnFirstFoundAdapter, false);

        if ( args.no5p3pAdapters ) {
            adapter = "$";
        } else if ( args.adapter3pOverride != null ) {
            adapter = args.adapter3pOverride;
        } else {
            if ( args.libraryDirection == SingleCellPipelineToolArgumentCollection.LibraryDirection.ThreePrime) {
                adapter = args.guide
                        ? "AGATCGGAAGAGCACACGTCTG" : "CCCATGTACTCTGCGTTGATACCACTGCTT";
            } else {
                adapter = args.guide
                        ? "CTGTCTCTTATACACATCT" : "AGATCGGAAGAGCACACGTCTG";
            }
        }
        adapter3pPattern = new AdapterUtils.AdapterPattern(adapter, args.adapterMinErrorRate, args.adapterMinOverlap, args.returnFirstFoundAdapter, true);

        if ( args.adapterMiddleOverride != null ) {
            adapter = args.adapterMiddleOverride;
        } else {
            if ( args.libraryDirection == SingleCellPipelineToolArgumentCollection.LibraryDirection.FivePrime ) {
                adapter = "^TTTCTTATATGGG";
            } else {
                adapter = args.guide
                        ? "GCTGTTTCCAGCTTAGCTCTTAAAC" : "XTTTTTTTTTTTTTTTTTTTTTTTTT";
            }
        }
        adapterMiddlePattern = new AdapterUtils.AdapterPattern(adapter, args.adapterMinErrorRate, args.adapterMinOverlap, args.returnFirstFoundAdapter, false);

        // other
        umiLengthOrig = (args.chemistry == SingleCellPipelineToolArgumentCollection.Chemistry.TenX_V3) ? 12 : 10;
        cbcUmiLength = args.umiLength + CBC_LENGTH;
        cbcUmiLengthOrig = umiLengthOrig + CBC_LENGTH;

        return super.instanceMainPostParseArgs();
    }

    private void logAdapters(final GATKRead read, FoundAdapters foundAdapters) {
        logAdapters(read, new AdapterUtils.FoundAdapter[] {foundAdapters.adapter5p, foundAdapters.adapterMiddle, foundAdapters.adapter3p});
    }

    private void logAdapters(final GATKRead read, final AdapterUtils.FoundAdapter[] foundAdapters) {

        // log read lines
        logger.info("N " + read.getName());
        logger.info("R " + read.getBasesString());

        // log adapters
        for ( final AdapterUtils.FoundAdapter fa : foundAdapters ) {
            final StringBuilder sb = new StringBuilder();
            sb.append(" ");
            if ( fa != null ) {
                for (int i = 0; i < fa.start; i++) {
                    sb.append(" ");
                }
                sb.append(">");
                for (int i = 0; i < fa.length; i++) {
                    sb.append(fa.adapter != null ? (char) fa.adapter.getPattern()[i] : '.');
                }
                sb.append("<");
            } else {
                sb.append("Not found");
            }
            logger.info(sb);
        }
    }

    @Override
    protected boolean shouldExitEarly(GATKRead read) {

        // limit number of input and output reads
        return ((args.maxInputReads != 0) && (stats.readsIn >= args.maxInputReads))
                || ((args.maxOutputReads != 0) && (stats.readsOut >= args.maxOutputReads));
    }
}
