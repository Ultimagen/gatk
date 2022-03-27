package org.broadinstitute.hellbender.tools.walkers.groundtruth;

import com.opencsv.CSVReader;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.FlowBasedProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.FlowBasedAlignmentArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.featuremapping.FlowFeatureMapper;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.FlowBasedAlignmentLikelihoodEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.LikelihoodEngineArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.haplotype.FlowBasedHaplotype;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.*;

import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.*;
import java.util.zip.GZIPOutputStream;

@CommandLineProgramProperties(
        summary = "Ground Truth Scorer",
        oneLineSummary = "Score reads against a reference/ground truth",
        programGroup = FlowBasedProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public class GroundTruthScorer extends ReadWalker {
    private static final Logger logger = LogManager.getLogger(GroundTruthScorer.class);
    public static final String OUTPUT_CSV_LONG_NAME = "output-csv";
    public static final String QUAL_REPORT_CSV_LONG_NAME = "qual-report-csv";
    public static final String HMER_REPORT_CSV_LONG_NAME = "hmer-report-csv";
    public static final String USE_SOFTCLIPPED_BASES_LONG_NAME = "use-softclipped-bases";
    public static final String GENOME_PRIOR_LONG_NAME = "genome-prior";
    public static final String FEATURES_FILE_LONG_NAME = "features-file";

    private static final int QUAL_VALUE_MAX = 60;
    private static final int HMER_VALUE_MAX = FlowBasedRead.MAX_CLASS;

    private static class ReportEntry {
        final int id;
        long count;
        double errorSum;
        double errorMin = Double.MAX_VALUE;
        double errorMax = Double.MIN_VALUE;

        ReportEntry(final int id) {
            this.id = id;
        }

        void add(final double error) {
            count++;
            errorSum += error;
            errorMin = Math.min(errorMin, error);
            errorMax = Math.max(errorMax, error);
        }

        String toCsvString() {
            if ( count == 0 ) {
                return id + ",0,0.0";
            } else {
                return id + String.format(",%d,%f", count, errorSum / count);
            }
        }

        static String csvHeading(final String idName) {
            return idName + ",count,value";
        }

        static ReportEntry[] newReport(final int size) {
            ReportEntry[]   report = new ReportEntry[size];
            for ( byte i = 0 ; i < report.length ; i++ ) {
                report[i] = new ReportEntry(i);
            }
            return report;
        }
    }


    @Argument(fullName = OUTPUT_CSV_LONG_NAME, doc="main CSV output file. supported file extensions: .csv, .csv.gz.")
    public GATKPath outputCsvPath = null;

    @Argument(fullName = QUAL_REPORT_CSV_LONG_NAME, doc="quality report output file. supported file extensions: .csv", optional = true)
    public GATKPath qualReportCsvPath = null;

    @Argument(fullName = HMER_REPORT_CSV_LONG_NAME, doc="hmer report output file. supported file extensions: .csv", optional = true)
    public GATKPath hmerReportCsvPath = null;

    @ArgumentCollection
    public LikelihoodEngineArgumentCollection likelihoodArgs = new LikelihoodEngineArgumentCollection();

    @ArgumentCollection
    public FlowBasedAlignmentArgumentCollection fbargs = new FlowBasedAlignmentArgumentCollection();

    @Argument(fullName = USE_SOFTCLIPPED_BASES_LONG_NAME, doc="", optional = true)
    public boolean useSoftclippedBases;

    @Argument(fullName = GENOME_PRIOR_LONG_NAME, doc="CSV input file containing genome-prior (one line per base with hmer frequencies).", optional = true)
    public GATKPath genomePriorPath;

    @Argument(fullName = FEATURES_FILE_LONG_NAME, doc="A VCF file containing features to be used as a use for filtering reads.", optional = true)
    public FeatureDataSource<VariantContext> features;

    // locals
    private FlowBasedAlignmentLikelihoodEngine likelihoodCalculationEngine;
    private PrintWriter                         outputCsv;
    private DecimalFormat                       doubleFormat = new DecimalFormat("0.0#####");
    private GenomePriorDB                       genomePriorDB;
    private ReportEntry[]                       qualReport;
    private ReportEntry[]                       hmerReport;

    // static/const
    static final private String[]       CSV_FIELD_ORDER = {
            "ReadName", "ReadKey", "ReadIsReversed", "ReadMQ", "ReadRQ", "GroundTruthKey", "ReadSequence", "Score", "ErrorProbability",
            "ReadKeyLength", "GroundTruthKeyLength", "CycleSkipStatus", "Cigar"
    };

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        // create likelihood engine
        final ReadLikelihoodCalculationEngine engine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(likelihoodArgs, false);
        if ( engine instanceof FlowBasedAlignmentLikelihoodEngine ) {
            likelihoodCalculationEngine = (FlowBasedAlignmentLikelihoodEngine)engine;
        } else {
            throw new GATKException("must use a flow based likelihood calculation engine");
        }

        // open genome prior if provided
        if ( genomePriorPath != null ) {
            try {
                genomePriorDB = new GenomePriorDB(genomePriorPath);
            } catch (IOException e) {
                throw new GATKException("failed to open genome-prior file: " + genomePriorPath);
            }
        }

        // open output, write header
        try {
            if (outputCsvPath.toPath().toString().endsWith(".gz")) {
                outputCsv = new PrintWriter(new GZIPOutputStream(outputCsvPath.getOutputStream()));
            } else {
                outputCsv = new PrintWriter(outputCsvPath.getOutputStream());
            }
        } catch (IOException e) {
            throw new GATKException("failed to open csv output: " + outputCsvPath, e);
        }
        emitCsvHeaders();

        // initialize reports
        if ( qualReportCsvPath != null ) {
            qualReport = ReportEntry.newReport(QUAL_VALUE_MAX + 1);
        }
        if ( hmerReportCsvPath != null ) {
            hmerReport = ReportEntry.newReport(HMER_VALUE_MAX + 1);
        }
    }

    @Override
    public void closeTool() {

        // close main output csv
        if ( outputCsv != null ) {
            outputCsv.close();
        }

        // write reports
        if ( qualReport != null ) {
            writeReport(qualReport, "qual", qualReportCsvPath);
        }
        if ( hmerReport != null ) {
            writeReport(hmerReport, "hmer", hmerReportCsvPath);
        }

        super.closeTool();
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        // handle read clipping
        final GATKRead clippedRead;
        if (ReadUtils.isSoftClipped(read) ) {
            if (useSoftclippedBases) {
                referenceContext.setWindow(read.getStart() - read.getUnclippedStart(), read.getUnclippedEnd() - read.getEnd());;
                clippedRead = ReadClipper.revertSoftClippedBases(read);
            } else {
                clippedRead = ReadClipper.hardClipSoftClippedBases(read);
            }
        } else {
            clippedRead = read;
        }

        // filter?
        if ( (features != null) && !filter(clippedRead, referenceContext) ) {
            return;
        }

        // create flow read/haplotype
        final FlowBasedReadUtils.ReadGroupInfo rgInfo = FlowBasedReadUtils.getReadGroupInfo(getHeaderForReads(), clippedRead);
        final FlowBasedRead flowRead = new FlowBasedRead(clippedRead, rgInfo.flowOrder, rgInfo.maxClass, fbargs);
        final Haplotype haplotype = new Haplotype(referenceContext.getBases(), true);
        final FlowBasedHaplotype flowHaplotype = new FlowBasedHaplotype(haplotype, rgInfo.flowOrder);

        // is this really needed?
        if ( !flowRead.isValid() ) {
           return;
        }

        // compute score
        final int         hapKeyLength = flowHaplotype.getKeyLength();
        final double      score = FlowFeatureMapper.computeLikelihoodLocal(flowRead, flowHaplotype, hapKeyLength, false);

        // compute error probability
        final double[]    errorProb = computeErrorProb(flowRead, genomePriorDB);

        // cycle skip
        final FlowBasedReadUtils.CycleSkipStatus cycleSkipStatus = FlowBasedReadUtils.getCycleSkipStatus(flowRead, referenceContext);

        // accumulate reports
        if ( hmerReport != null ) {
            addToHmerReport(flowRead.getKey(), errorProb);
        }
        if ( cycleSkipStatus == FlowBasedReadUtils.CycleSkipStatus.NS &&  qualReport != null ) {
            addToQualReport(flowRead, referenceContext, errorProb);
        }

        // emit
        try {
            emit(flowRead, flowHaplotype, score, errorProb, read, referenceContext, cycleSkipStatus);
        } catch (IOException e) {
            throw new GATKException("failed to write output record", e);
        }
    }

    private boolean filter(final GATKRead read, final ReferenceContext referenceContext) {

        // loop on features contained in the read, check that they are in agreement with read data
        Iterator<VariantContext>    iter = features.query(new SimpleInterval(read));
        byte[]                      ref = referenceContext.getBases();
        while ( iter.hasNext() ) {
            final VariantContext vc = iter.next();
            for ( int refCoord = vc.getStart() ; refCoord <= vc.getEnd() ; refCoord++ ) {

                // get byte from read
                Optional<Byte>      readByte = ReadUtils.getReadBaseAtReferenceCoordinate(read, refCoord);
                if ( !readByte.isPresent() ) {
                    return false;
                }

                // get byte from reference
                byte                refByte = ref[refCoord - referenceContext.getWindow().getStart()];

                // compare
                if ( refByte != readByte.get() ) {
                    return false;
                }
            }
        }

        // if here, no interference detected
        return true;
    }

    /*
     * compute error probability vector for a read
     *
     * The vector has one element for each flow key, representing the probability complementing the call-probability to 1
     * This is further complicated by the optional presence of a genome-prior database, which provides factoring for
     * each hmer length (on a base basis)
     */
    private double[] computeErrorProb(final FlowBasedRead flowRead, final GenomePriorDB genomePriorDB) {

        final int[] key = flowRead.getKey();
        final byte[] flowOrder = flowRead.getFlowOrderArray();
        final double[] probCol = new double[FlowBasedRead.MAX_CLASS + 1];
        double[] result = new double[key.length];

        for ( int i = 0 ; i < key.length ; i++ ) {

            // step 1 - extract column & sum
            double  sum = 0;
            for ( int j = 0 ; j < probCol.length ; j++ ) {
                sum += (probCol[j] = flowRead.getProb(i, j));
            }

            // step 2 - normalize column
            for ( int j = 0 ; j < probCol.length ; j++ ) {
                probCol[j] /= sum;
            }

            // step 3 - scale according to prior genome?
            if ( genomePriorDB != null ) {

                long[] prior = (genomePriorDB != null) ? genomePriorDB.getPriorForBase(flowOrder[i]) : null;
                sum = 0;
                for ( int j = 0 ; j < probCol.length ; j++ ) {
                    sum += (probCol[j] *= prior[j]);
                }

                // assign normalized result
                result[i] = 1 - (probCol[Math.min(key[i], flowRead.getMaxHmer())] / sum);
            } else {

                // assign normalized result
                result[i] = 1 - probCol[Math.min(key[i], flowRead.getMaxHmer())];
            }
        }

        return result;
    }

    private void emitCsvHeaders() {

        outputCsv.println(StringUtils.join(CSV_FIELD_ORDER, ","));
    }

    private void emit(final FlowBasedRead flowRead, final FlowBasedHaplotype refHaplotype, final double score, final double[] errorProb,
                      GATKRead read, ReferenceContext referenceContext,
                      FlowBasedReadUtils.CycleSkipStatus cycleSkipStatus) throws IOException {

        // build line columns
        final Map<String,Object> cols = new LinkedHashMap<>();

        // read info
        cols.put("ReadName", flowRead.getName());
        cols.put("ReadIsReversed", flowRead.isReverseStrand() ? 1 : 0);
        cols.put("ReadMQ", flowRead.getMappingQuality());
        cols.put("ReadRQ", flowRead.getAttributeAsFloat("rq"));
        cols.put("CycleSkipStatus", cycleSkipStatus);
        cols.put("Cigar", read.getCigar().toString());


        // keys, seq, etc
        cols.put("ReadKey", "\"" + StringUtils.join(flowRead.getKey(), ',') + "\"");
        cols.put("GroundTruthKey", "\"" + StringUtils.join(refHaplotype.getKey(), ',') + "\"");
        cols.put("ReadSequence", flowRead.getBasesString());
        cols.put("ReadKeyLength", flowRead.getKeyLength());
        cols.put("GroundTruthKeyLength", refHaplotype.getKeyLength());

        // scores
        cols.put("Score", score);
        cols.put("ErrorProbability", "\"" + StringUtils.join(
                Arrays.stream(errorProb).mapToObj(v -> doubleFormat.format(v)).toArray(),
                ',') + "\"");

        // construct line
        StringBuilder       sb = new StringBuilder();
        int                 colIndex = 0;
        for ( String field : CSV_FIELD_ORDER ) {
            if ( colIndex++ > 0 ) {
                sb.append(',');
            }
            if ( !cols.containsKey(field) ) {
                throw new GATKException("column missing from csv line: " + field);
            }
            sb.append(cols.get(field));
            cols.remove(field);
        }
        if ( cols.size() > 0 ) {
            throw new GATKException("invalid columns on csv line: " + cols.keySet());
        }

        // output line
        outputCsv.println(sb);
    }

    static private class GenomePriorDB {

        final private Map<Byte, long[]>      db = new LinkedHashMap<>();

        GenomePriorDB(GATKPath path) throws IOException {

            final CSVReader     csvReader = new CSVReader(new InputStreamReader(path.getInputStream()));
            String[]            line;
            while ( (line = csvReader.readNext()) != null ) {
                long[]          prior = new long[FlowBasedRead.MAX_CLASS + 1];
                Byte            base = line[0].getBytes()[0];
                for ( int i = 0 ; i < prior.length ; i++ ) {
                    prior[i] = Long.parseLong(line[i+1]);
                }
                db.put(base, prior);
            }
        }

        long[]   getPriorForBase(byte base) {
            return db.get(base);
        }
    }

    private void addToHmerReport(final int[] key, final double[] errorProb) {

        for ( int i = 0 ; i < key.length ; i++ ) {
            if ( key[i] < hmerReport.length ) {
                hmerReport[key[i]].add(errorProb[i]);
            }
        }
    }

    private void writeReport(final ReportEntry[] report, final String idName, final GATKPath outputCsvPath) {

        PrintWriter     pw = new PrintWriter(outputCsvPath.getOutputStream());
        pw.println(ReportEntry.csvHeading(idName));
        for ( int i = 0 ; i < report.length ; i++ ) {
            pw.println(report[i].toCsvString());
        }
        pw.close();
    }

    private void addToQualReport(FlowBasedRead flowRead, ReferenceContext referenceContext, final double[] errorProb) {

        // convert reference to key space
        final Haplotype             haplotype = new Haplotype(referenceContext.getBases(), true);
        final FlowBasedHaplotype    flowHaplotype = new FlowBasedHaplotype(haplotype, flowRead.getFlowOrder());

        // access keys
        final int[]                 readKey = flowRead.getKey();
        final int[]                 hapKey = flowHaplotype.getKey();

        // loop on key positions
        for ( int flow = 0 ; flow < Math.min(readKey.length, hapKey.length) ; flow++ ) {

            // determine quality
            final double        prob = errorProb[flow];
            final int           qual = (int)(-10 * Math.log10(prob));

            // determine if matches reference
            final boolean       diff = readKey[flow] != hapKey[flow];

            // accumulate
            if ( qual < qualReport.length ) {
                qualReport[qual].add(diff ? 1.0 : 0.0);
            }
        }
    }


}
