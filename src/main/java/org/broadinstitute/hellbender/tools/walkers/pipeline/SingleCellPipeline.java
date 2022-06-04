package org.broadinstitute.hellbender.tools.walkers.pipeline;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.util.Tuple;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.Arrays;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

public class SingleCellPipeline {

    private static final Logger logger = LogManager.getLogger(SingleCellPipeline.class);
    public static final int CBC_LENGTH = 16;
    public static final int CDNA_MAX_LENGTH = 315;
    public static final String REPORT_JSON_FILEnAME_SUFFIX = "_report.json";
    public static final char CBC_UMI_MASK_BASE = 'T';
    public static final char CBC_UMI_MASK_QUAL = 40;

    // public arguments
    final public SingleCellPipelineToolArgumentCollection args;

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
    private CbcWhitelist cbcWhitelist;

    private int umiLengthOrig;
    private int cbcUmiLength;
    private int cbcUmiLengthOrig;

    private BlockingQueue<Tuple<FastqRecord, FastqRecord>> outputQueue;
    private Thread outputThread;


    public String lastReason;

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

    public SingleCellPipeline(final SingleCellPipelineToolArgumentCollection args) {
        this.args = args;
    }

    public void process(final String readName, final boolean isReverseStrand, byte[] bases, byte[] quals, final AttributeProvider attributeProvider) {

        final boolean debug = readName.equals(args.debugTrimmingFor);
        stats.readsIn++;
        stats.bpIn += bases.length;

        // rsq threshold?
        if ( args.rsqThreshold > 0 && attributeProvider.hasAttribute("RQ")
            && attributeProvider.getAttributeAsInteger("RQ") < args.rsqThreshold ) {
            stats.rqDropped++;
            lastReason = "rsq threshold";
            return;
        }

        // reverse read?
        if ( isReverseStrand ) {
            bases = BaseUtils.simpleReverseComplement(bases);
            quals = Arrays.copyOf(quals, quals.length);
            ArrayUtils.reverse(quals);
        }

        // access read
        final int basesTrimmedLength;
        basesTrimmedLength = findTrimmingPoint(readName, bases, quals, 0, quals.length, args.qualityCutoff);
        stats.bpCutoff += (bases.length - basesTrimmedLength);
        if (basesTrimmedLength < (cbcUmiLength + args.minCdnaLength)) {
            stats.trimmedTooShort++;
            lastReason = "too short after trimming";
            return;
        }

        // find adapters
        FoundAdapters foundAdapters = findAdapters(readName, bases, quals, basesTrimmedLength, debug);

        if ( foundAdapters != null ) {

            // find and limit read1 (cbc_umi)
            int read1Start = foundAdapters.adapter5p.start + foundAdapters.adapter5p.length;
            int read1Length = Math.min(foundAdapters.adapter3p.start - read1Start, cbcUmiLengthOrig);
            boolean read1Valid = read1Length == cbcUmiLengthOrig;
            stats.read1TooShortDropped += (read1Valid ? 0 : 1);
            if ( !read1Valid ) {
                lastReason = "read1 too short";
            }

            // find and limit read1 (cdna)
            int read2Start = foundAdapters.adapterMiddle.start + foundAdapters.adapterMiddle.length;
            int read2Length = Math.min(foundAdapters.adapter3p.start - read2Start, CDNA_MAX_LENGTH);
            boolean read2Valid = read2Length >= args.minCdnaLength;
            if ( !read2Valid ) {
                lastReason = "read2 too short: " + read2Length  + " vs " + args.minCdnaLength;
            }
            if ( read2Valid ) {
                if (args.cdnaFirstBasesToClip != 0) {
                    read2Start += args.cdnaFirstBasesToClip;
                    read2Length -= args.cdnaFirstBasesToClip;
                }
                if (args.cdnaTrimmingLength != 0) {
                    read2Length = Math.min(read2Length, args.cdnaTrimmingLength);
                }
                if (read2Length <= 0) {
                    read2Valid = false;
                    if ( !read2Valid ) {
                        lastReason = "read2 too short after cdnaTrimmingLength";
                    }
                }
            }
            stats.read2TooShortDropped += (read2Valid ? 0 : 1);
            if ( debug ) {
                logger.info("read2start: " + read2Start + ", read2Length: " + read2Length + ", read2Valid: " + read2Valid + ", reason:" + lastReason);
            }

            // filter for umi quality
            if ( read1Valid ) {
                if (args.umiQualityThreshold > 0) {
                    if (anyQualBelowThreshold(quals, read1Start + read1Length - umiLengthOrig, umiLengthOrig, args.umiQualityThreshold)) {
                        stats.umiQualityDropped++;
                        read1Valid = false;
                        lastReason = "read1Valid turned false because umiQualityThreshold";
                    }
                }
            }

            // match to whitelist or correct
            byte[] read1bases = bases;
            if ( read1Valid ) {
                if (cbcWhitelist != null) {
                    stats.processedCbcWhitelist++;
                    CbcWhitelist.Result result = cbcWhitelist.matchOrCorrect(bases, read1Start, CBC_LENGTH);
                    switch (result.type) {
                        case NOT_FOUND:
                            stats.notFoundCbcWhitelist++;
                            break;
                        case EXACT:
                            stats.matchedCbcWhitelist++;
                            break;
                        case AMBIGUOUS:
                            stats.ambiguousCbcWhitelist++;
                            break;
                        case DEL_CORRECTED:
                            stats.delCorrectedCbcWhitelist++;
                            read1bases = new byte[read1Length + 1];
                            Utils.validate(result.correction.length() == CBC_LENGTH, "DEL correction same as CBC length");
                            System.arraycopy(result.correction.getBytes(), 0, read1bases, 0, CBC_LENGTH);
                            System.arraycopy(bases, read1Start + CBC_LENGTH - 1, read1bases, CBC_LENGTH, read1Length - CBC_LENGTH + 1);
                            read1Start = 0;
                            read1Length++;
                            break;
                        case INS_CORRECTED:
                            stats.insCorrectedCbcWhitelist++;
                            read1bases = new byte[read1Length - 1];
                            Utils.validate(result.correction.length() == CBC_LENGTH, "INS correction same as CBC length");
                            System.arraycopy(result.correction.getBytes(), 0, read1bases, 0, CBC_LENGTH);
                            System.arraycopy(bases, read1Start + CBC_LENGTH + 1, read1bases, CBC_LENGTH, read1Length - CBC_LENGTH - 1);
                            read1Start = 0;
                            read1Length--;
                            break;
                        case SNP_CORRECTED:
                            stats.snpCorrectedCbcWhitelist++;
                            read1bases = new byte[read1Length];
                            Utils.validate(result.correction.length() == CBC_LENGTH, "INS correction same as CBC length");
                            System.arraycopy(result.correction.getBytes(), 0, read1bases, 0, CBC_LENGTH);
                            System.arraycopy(bases, read1Start + CBC_LENGTH, read1bases, CBC_LENGTH, read1Length - CBC_LENGTH - 1);
                            read1Start = 0;
                            break;
                    }
                }
            }

            // reads valid?
            if ( read1Valid && read2Valid ) {

                stats.readsOut++;
                stats.bpOut += (read1Length + read2Length);

                if (args.logAdapters == SingleCellPipelineToolArgumentCollection.LogAdapters.Output) {
                    logAdapters(readName, bases, quals, foundAdapters);
                }

                if ( !args.noOutput ) {
                    final FastqRecord rec1 = makeFastQRecord(readName, read1bases, quals, read1Start, read1Start + read1Length, false, true);
                    final FastqRecord rec2 = makeFastQRecord(readName, bases, quals, read2Start, read2Start + read2Length, args.reverseComplementRead2, false);

                    if ( args.multiproc ) {
                        outputQueue.offer(new Tuple<>(rec1, rec2));
                    } else {
                        read1Writer.write(rec1);
                        read2Writer.write(rec2);
                    }
                }
            } else {
                if ( lastReason == null ) {
                    lastReason = "no output (general)";
                }
            }
        }
    }

    private FoundAdapters findAdapters(String readName, byte[] bases, byte[] quals, int length, boolean debug) {

        FoundAdapters result = null;

        // look for the adapters
        final AdapterUtils.FoundAdapter adapter5p = AdapterUtils.findAdapter(bases, adapter5pPattern, 0, length);
        final AdapterUtils.FoundAdapter adapterMiddle = AdapterUtils.findAdapter(bases, adapterMiddlePattern, 0, length);
        final AdapterUtils.FoundAdapter adapter3p;
        if ( adapterMiddle != null ) {
            adapter3p = AdapterUtils.findAdapter(bases, adapter3pPattern,
                    adapterMiddle.start + adapterMiddle.length, length - adapterMiddle.start - adapterMiddle.length);
        } else {
            adapter3p = AdapterUtils.findAdapter(bases, adapter3pPattern, 0, length);
        }

        // stats
        stats.adapter5p += (adapter5p != null ? 1 : 0);
        stats.adapterMiddle += (adapterMiddle != null ? 1 : 0);
        stats.adapter3p += (adapter3p != null ? 1 : 0);

        // generate result. if adapter3p is not found, fallback on the end of the read
        if ( adapter5p != null && adapterMiddle != null ) {
            result = new FoundAdapters(adapter5p, adapterMiddle,
                    adapter3p != null ? adapter3p : new AdapterUtils.FoundAdapter(length, 0));
        }

        // log adapters?
        if ( debug || args.logAdapters == SingleCellPipelineToolArgumentCollection.LogAdapters.Input ) {
            logAdapters(readName, bases, quals,
                    new AdapterUtils.FoundAdapter[] {adapter5p, adapterMiddle, adapter3p});
        }

        return result;

    }

    /*
     * implementation of a 3' cutoff similar to cutadapt
     *
     * see https://cutadapt.readthedocs.io/en/stable/algorithms.html for details
     */
    private int findTrimmingPoint(String readName, final byte[] bases, final byte[] quals, final int start, final int end, int qualityCutoff) {

        final boolean debug = readName.equals(args.debugTrimmingFor);
        final boolean nextGen = false;

        // debugging trimming for this read?
        if ( debug ) {
            logger.info("Trimming debug for: " + readName + " quals.length: " + quals.length
                                                            + ", start: " + start + ", end: " + end);
        }

        // sanity
        Utils.validate(end >= start, "end can't be before start");
        if ( end == start ) {
            if ( debug ) {
                logger.info("end same as start");
            }
            return start;
        }

        // create partial sums of quality values reduces by threshold. Stop when positive. stop if greater than zero
        int q = (nextGen && bases[end-1] == 'G') ? (qualityCutoff - 1) : quals[end - 1];
        int sum = q - qualityCutoff;
        int minIndex = end - 1;
        int minValue = sum;
        if ( debug ) {
            logger.info("" + (end - 1) + ": " + quals[end-1] + " sum: " + sum);
        }
        for ( int i = end - 2 ; i >= start ; i-- ) {
            q = (nextGen && bases[i] == 'G') ? (qualityCutoff - 1) : quals[i];
            sum += (q - qualityCutoff);
            if ( debug ) {
                logger.info("" + i + ": " + quals[i] + " sum: " + sum);
            }
            if ( sum < minValue ) {
                minIndex = i;
                minValue = sum;
            }
            if ( sum > 0 )
                break;
        }

        // trim?
        if ( minValue < 0) {
            if ( debug ) {
                logger.info("trimmed at minIndex of: " + minIndex + ", length will be: " + (minIndex - start));
            }
            return Math.max(0, minIndex);
        } else {
            if ( debug ) {
                logger.info("trimmed at end (not trimmed): " + end + ", length will be: " + (minIndex - start));
            }
            return end;
        }

    }

    private boolean anyQualBelowThreshold(final byte[] quals, final int start, final int length, final int threshold) {
        for ( int i = 0 ; i < length ; i++ ) {
            if ( quals[start + i] < threshold ) {
                return true;
            }
        }

        return false;
    }

    private FastqRecord makeFastQRecord(String readName, byte[] readBases, byte[] readQuals, int startOfs, int endOfs, boolean rc, boolean maskLast) {

        final int length = endOfs - startOfs;
        Utils.validIndex(startOfs, readBases.length, "startOfs: out of range");
        Utils.validate(startOfs + length <= readBases.length, "startOfs + length: out of range");
        byte[] bases = Arrays.copyOfRange(readBases, startOfs, startOfs + length);
        final byte[] quals = Arrays.copyOfRange(readQuals, startOfs, startOfs + length);
        if ( maskLast && args.cbcUmiMaskLastBytes > 0 ) {
            for ( int i = 0 ; i < args.cbcUmiMaskLastBytes ; i++ ) {
                bases[bases.length - i - 1] = CBC_UMI_MASK_BASE;
                quals[quals.length - i - 1] = CBC_UMI_MASK_QUAL;
            }
        }
        if ( rc ) {
            ArrayUtils.reverse(quals);
            bases = BaseUtils.simpleReverseComplement(bases);
        }


        return new FastqRecord(readName, bases, null, quals);
    }

    public void onStart() {

        // log adapters & other key parameters
        logger.info("args.baseFilename: " + args.baseFilename);
        logger.info("args.guide: " + args.guide);
        logger.info("args.libraryDirection: " + args.libraryDirection);
        logger.info("args.chemistry: " + args.chemistry);

        logger.info("args.illumina: " + args.illumina);

        logger.info("args.umiLength: " + args.umiLength);
        logger.info("umiLengthOrig: " + umiLengthOrig);
        logger.info("cbcUmiLength: " + cbcUmiLength);
        logger.info("cbcUmiLengthOrig: " + cbcUmiLengthOrig);
        logger.info("args.cbcUmiMaskLastBytes: " + args.cbcUmiMaskLastBytes);

        logger.info("args.minCdnaLength: " + args.minCdnaLength);
        logger.info("args.cdnaFirstBasesToClip: " + args.cdnaFirstBasesToClip);
        logger.info("args.cdnaTrimmingLength: " + args.cdnaTrimmingLength);

        logger.info("args.rsqThreshold: " + args.rsqThreshold);
        logger.info("args.illuminaRead1List: " + args.illuminaRead1List);
        logger.info("args.illuminaRead2List: " + args.illuminaRead2List);

        logger.info("args.adapterMinErrorRate: " + args.adapterMinErrorRate);
        logger.info("args.adapterMinOverlap: " + args.adapterMinOverlap);
        logger.info("args.qualityCutoff: " + args.qualityCutoff);
        logger.info("args.umiQualityThreshold: " + args.umiQualityThreshold);
        logger.info("args.reverseComplementRead2: " + args.reverseComplementRead2);

        logger.info("args.no5p3pAdapters: " + args.no5p3pAdapters);
        logger.info("adapter5p: " + adapter5pPattern.getDescription());
        logger.info("adapter3p: " + adapter3pPattern.getDescription());
        logger.info("adapterMiddle: " + adapterMiddlePattern.getDescription());

        logger.info("args.cbcWhitelistPath: " + args.cbcWhitelistPath);

        logger.info("args.noOutput: " + args.noOutput);
        logger.info("args.logAdapters: " + args.logAdapters);
        logger.info("args.maxInputReads: " + args.maxInputReads);
        logger.info("args.maxOutputReads: " + args.maxOutputReads);

        // open writers
        File f;
        fastqWriterFactory.setCreateMd5(false);
        read1Writer = fastqWriterFactory.newWriter(f = buildFastQReedsOutputFile(1));
        logger.info("read1 output: " + f);
        read2Writer = fastqWriterFactory.newWriter(f = buildFastQReedsOutputFile(2));
        logger.info("read2 output: " + f);

        // read whitelist
        if ( args.cbcWhitelistPath != null ) {
            cbcWhitelist = new CbcWhitelist(args.cbcWhitelistPath, args.cbcWhitelistSupportsSnp);
        }

        // create queue
        if ( args.multiproc ) {
            outputQueue = new LinkedBlockingQueue<>(args.queueCapacity);
            (outputThread = new Thread(() -> {
                try {
                    while (true) {
                        final Tuple<FastqRecord,FastqRecord> tuple = outputQueue.take();
                        if (tuple.a == null) {
                            break;
                        } else {
                            read1Writer.write(tuple.a);
                            read2Writer.write(tuple.b);
                        }
                    }
                } catch (InterruptedException e) {
                    logger.warn("", e);
                }
            })).start();
        }



    }

    private File buildFastQReedsOutputFile(int i) {
        return new File(args.baseFilename + "_read" + i + ".fastq" + (args.compressedOutput ? ".gz" : ""));
    }

    public void onClose() {

        // wait for writer
        if ( args.multiproc ) {
            outputQueue.offer(new Tuple<>(null, null));
            try {
                outputThread.join();
            } catch (InterruptedException e) {
                logger.warn("", e);
            }
            ;
        }

        // close writers
        if ( read1Writer != null ) {
            read1Writer.close();
        }
        if ( read2Writer != null ) {
            read2Writer.close();
        }

        // write stats
        stats.writeJson(new File(args.baseFilename + REPORT_JSON_FILEnAME_SUFFIX), args);
    }

    // perform additional argument verification and adjustments, assign defaults
    public void onArgsReady() {

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
        adapter5pPattern = new AdapterUtils.AdapterPattern(adapter, args.adapterMinErrorRate, args.adapterMinOverlap, true, false);

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
        adapter3pPattern = new AdapterUtils.AdapterPattern(adapter, args.adapterMinErrorRate, args.adapterMinOverlap, true, true);

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
        adapterMiddlePattern = new AdapterUtils.AdapterPattern(adapter, args.adapterMinErrorRate, args.adapterMinOverlap, false, false);

        // other
        umiLengthOrig = (args.chemistry == SingleCellPipelineToolArgumentCollection.Chemistry.V3) ? 12 : 10;
        cbcUmiLength = args.umiLength + CBC_LENGTH;
        cbcUmiLengthOrig = umiLengthOrig + CBC_LENGTH;
    }

    private void logAdapters(final String readName, final byte[] bases, final byte[] quals, final FoundAdapters foundAdapters) {
        logAdapters(readName, bases, quals, new AdapterUtils.FoundAdapter[] {foundAdapters.adapter5p, foundAdapters.adapterMiddle, foundAdapters.adapter3p});
    }

    private void logAdapters(final String readName, final byte[] bases, final byte[] quals, final AdapterUtils.FoundAdapter[] foundAdapters) {

        // log read lines
        logger.info("N " + readName);
        logger.info("Q " + quals2String(quals));
        logger.info("R " + new String(bases));

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
                    sb.append((fa.adapter != null) && (i < fa.adapter.getPattern().length) ? (char) fa.adapter.getPattern()[i] : '.');
                }
                sb.append("<");
            } else {
                sb.append("Not found");
            }
            logger.info(sb);
        }
    }

    private String quals2String(byte[] quals) {
        StringBuilder   sb = new StringBuilder();
        for ( byte q : quals ) {
            sb.append((char)(q + '!'));
        }
        return sb.toString();
    }

    public boolean shouldExitEarly() {

        // limit number of input and output reads
        return ((args.maxInputReads != 0) && (stats.readsIn >= args.maxInputReads))
                || ((args.maxOutputReads != 0) && (stats.readsOut >= args.maxOutputReads));
    }

}
