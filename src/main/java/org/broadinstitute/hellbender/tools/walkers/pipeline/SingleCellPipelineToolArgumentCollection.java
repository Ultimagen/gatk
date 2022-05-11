package org.broadinstitute.hellbender.tools.walkers.pipeline;

import org.apache.commons.collections.CollectionUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.GATKException;

import java.io.Serializable;
import java.util.List;

public class SingleCellPipelineToolArgumentCollection implements Serializable {
    private static final long serialVersionUID = 0;
    public static final String LONG_NAME_BASE_FILENAME = "base-filename";
    public static final String LONG_NAME_UMI_LENGTH = "umi-length";
    public static final String LONG_NAME_GUIDE = "guide";
    public static final String LONG_NAME_LIBRARY_DIRECTION = "library-direction";
    public static final String LONG_NAME_CHEMISTRY = "chemistry";
    public static final String LONG_NAME_ADAPTER_5_P_OVERRIDE = "adapter-5p-override";
    public static final String LONG_NAME_ADAPTER_3_P_OVERRIDE = "adapter-3p-override";
    public static final String LONG_NAME_ADAPTER_MIDDLE_OVERRIDE = "adapter-middle-override";
    public static final String LONG_NAME_ILLUMINA = "illumina";
    public static final String LONG_NAME_ILLUMINA_READ_1_LIST = "illumina-read1-list";
    public static final String LONG_NAME_ILLUMINA_READ_2_LIST = "illumina-read2-list";
    public static final String LONG_NAME_ADAPTER_MIN_ERROR_RATE = "adapter-min-error-rate";
    public static final String LONG_NAME_ADAPTER_MIN_OVERLAP = "adapter-min-overlap";
    public static final String LONG_NAME_NO_5_P_3_P_ADAPTERS = "no-5p-3p-adapters";
    public static final String LONG_NAME_REVERSE_COMPLEMENT_READ_2 = "reverse-complement-read2";
    public static final String LONG_NAME_NO_OUTPUT = "no-output";
    public static final String LONG_NAME_RETURN_FIRST_FOUND_ADAPTER = "return-first-found-adapter";
    public static final String LONG_NAME_QUALITY_CUTOFF = "quality-cutoff";
    public static final String FULL_NAME_LOG_ADAPTERS = "log-adapters";
    public static final String FULL_NAME_MAX_OUTPUT_READS = "max-output-reads";
    public static final String FULL_NAME_MAX_INPUT_READS = "max-input-reads";
    public static final String LONG_NAME_MIN_CDNA_LENGTH = "min-cdna-length";
    public static final String LONG_NAME_CDNA_FIRST_BASES_TO_CLIP = "cdna-first-bases-to-clip";
    public static final String LONG_NAME_CDNA_TRIMMING_LENGTH = "cdna-trimming-length";
    public static final String LONG_NAME_RSQ_THRESHOLD = "rsq-threshold";
    public static final String FULL_NAME_CBC_UMI_MASK_LAST_BYTES = "cbc_umi_mask_last_bytes";

    @Argument(fullName = LONG_NAME_BASE_FILENAME, doc = "output files base name (prefix)")
    public String baseFilename;

    @Argument(fullName = LONG_NAME_GUIDE, doc = "is guide (true) or hash (false)?", optional = true)
    public boolean guide;

    @Argument(fullName = LONG_NAME_LIBRARY_DIRECTION, doc = "input library direction", optional = true)
    public LibraryDirection libraryDirection = LibraryDirection.ThreePrime;

    @Argument(fullName = LONG_NAME_CHEMISTRY, doc = "input chemistry", optional = true)
    public Chemistry chemistry = Chemistry.V3;

    @Argument(fullName = LONG_NAME_UMI_LENGTH, doc = "umi length", optional = true)
    public int umiLength = 10;

    @Argument(fullName = LONG_NAME_MIN_CDNA_LENGTH, doc = "minimal cdna length", optional = true)
    public int minCdnaLength = 20;

    @Argument(fullName = LONG_NAME_CDNA_FIRST_BASES_TO_CLIP, doc = "number of bases to clip from start of cdna", optional = true)
    public int cdnaFirstBasesToClip = 0;

    @Argument(fullName = LONG_NAME_CDNA_TRIMMING_LENGTH, doc = "length to trim cdna to", optional = true)
    public int cdnaTrimmingLength = 0;

    @Argument(fullName = FULL_NAME_CBC_UMI_MASK_LAST_BYTES, doc = "numbr of bytes to mask at the end of the cbc_umi", optional = true)
    public int cbcUmiMaskLastBytes;

    @Argument(fullName = LONG_NAME_RSQ_THRESHOLD, doc = "RQ attribute value to filter on", optional = true)
    public int rsqThreshold;

    @Argument(fullName = LONG_NAME_ADAPTER_5_P_OVERRIDE, doc = "override for the 5p adapter", optional = true)
    public String adapter5pOverride;

    @Argument(fullName = LONG_NAME_ADAPTER_3_P_OVERRIDE, doc = "override for the 3p adapter", optional = true)
    public String adapter3pOverride;

    @Argument(fullName = LONG_NAME_ADAPTER_MIDDLE_OVERRIDE, doc = "override for the middle adapter", optional = true)
    public String adapterMiddleOverride;

    @Argument(fullName = LONG_NAME_NO_5_P_3_P_ADAPTERS, doc = "input does not contain 5p/3p adapters (already trimmed)", optional = true)
    public boolean no5p3pAdapters;

    @Argument(fullName = LONG_NAME_REVERSE_COMPLEMENT_READ_2, doc = "read2 needs reverse complementing", optional = true)
    public boolean reverseComplementRead2 = true;

    @Argument(fullName = LONG_NAME_ILLUMINA, doc = "illumina mode (currently not supported)", optional = true)
    public boolean illumina;

    @Argument(fullName = LONG_NAME_ILLUMINA_READ_1_LIST, doc = "read1 files for illumina mode", optional = true)
    public List<GATKPath> illuminaRead1List;

    @Argument(fullName = LONG_NAME_ILLUMINA_READ_2_LIST, doc = "read2 filtes for illumina mode", optional = true)
    public List<GATKPath> illuminaRead2List;

    @Argument(fullName = LONG_NAME_ADAPTER_MIN_ERROR_RATE, doc = "error rate threshold for adapters", optional = true)
    public double adapterMinErrorRate = 0.2;

    @Argument(fullName = LONG_NAME_ADAPTER_MIN_OVERLAP, doc = "minimal overlap for adapters", optional = true)
    public int adapterMinOverlap = 10;

    @Argument(fullName = LONG_NAME_QUALITY_CUTOFF, doc = "quality cutoff value", optional = true)
    public int qualityCutoff = 30;

    @Argument(fullName = LONG_NAME_NO_OUTPUT, doc = "generate no output files (read1/read2)", optional = true)
    public boolean noOutput;

    @Argument(fullName = LONG_NAME_RETURN_FIRST_FOUND_ADAPTER, doc = "use first adapter found (rather than best)", optional = true)
    public boolean returnFirstFoundAdapter = true;

    @Argument(fullName = FULL_NAME_LOG_ADAPTERS, doc = "log adapters found on input or output reads", optional = true)
    public LogAdapters logAdapters = LogAdapters.None;

    @Argument(fullName = FULL_NAME_MAX_INPUT_READS, doc = "limit number of input reads processed", optional = true)
    public int maxInputReads;

    @Argument(fullName = FULL_NAME_MAX_OUTPUT_READS, doc = "limit numner of output reads generated", optional = true)
    public int maxOutputReads;

    enum LibraryDirection {
        ThreePrime,
        FivePrime
    }

    enum Chemistry {
        V2,
        V3
    }

    enum LogAdapters {
        None,
        Input,
        Output
    }

    protected void validate() {

        // illumina read lists apply to illumina mode only
        if ( !illumina ) {
            if (!CollectionUtils.isEmpty(illuminaRead1List)) {
                throw new IllegalArgumentException("--" + LONG_NAME_ILLUMINA_READ_1_LIST + " can only be used in illumina mode");
            }
            if (!CollectionUtils.isEmpty(illuminaRead2List)) {
                throw new IllegalArgumentException("--" + LONG_NAME_ILLUMINA_READ_2_LIST + " can only be used in illumina mode");
            }
        } else {
            throw new GATKException("illumina not supported");
        }
    }
}
