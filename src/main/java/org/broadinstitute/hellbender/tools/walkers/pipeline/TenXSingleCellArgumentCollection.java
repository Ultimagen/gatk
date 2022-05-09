package org.broadinstitute.hellbender.tools.walkers.pipeline;

import org.apache.commons.collections.CollectionUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.GATKPath;

import java.io.Serializable;
import java.util.List;

public class TenXSingleCellArgumentCollection implements Serializable {
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
    public static final String FULL_NANE_FASTQ_ASYNC_IO = "fastq-async-io";
    public static final String LONG_NAME_BYPASS_MODE = "bypass-mode";
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

    @Argument(fullName = LONG_NAME_BASE_FILENAME)
    public String baseFilename;

    @Argument(fullName = LONG_NAME_UMI_LENGTH, optional = true)
    public int umiLength = 10;

    @Argument(fullName = LONG_NAME_MIN_CDNA_LENGTH, optional = true)
    public int minCdnaLength = 20;

    @Argument(fullName = LONG_NAME_CDNA_FIRST_BASES_TO_CLIP, optional = true)
    public int cdnaFirstBasesToClip = 0;

    @Argument(fullName = LONG_NAME_CDNA_TRIMMING_LENGTH, optional = true)
    public int cdnaTrimmingLength = 0;

    @Argument(fullName = FULL_NAME_CBC_UMI_MASK_LAST_BYTES)
    public int cbcUmiMaskLastBytes;

    @Argument(fullName = LONG_NAME_RSQ_THRESHOLD)
    public int rsqThreshold;



    @Argument(fullName = LONG_NAME_GUIDE, doc = "=is guide (true) or hash(false)", optional = true)
    public boolean guide;

    @Argument(fullName = LONG_NAME_LIBRARY_DIRECTION, optional = true)
    public LibraryDirection libraryDirection = LibraryDirection.ThreePrime;

    @Argument(fullName = LONG_NAME_CHEMISTRY, optional = true)
    public Chemistry chemistry = Chemistry.TenX_V3;

    @Argument(fullName = LONG_NAME_ADAPTER_5_P_OVERRIDE, optional = true)
    public String adapter5pOverride;

    @Argument(fullName = LONG_NAME_ADAPTER_3_P_OVERRIDE, optional = true)
    public String adapter3pOverride;

    @Argument(fullName = LONG_NAME_ADAPTER_MIDDLE_OVERRIDE, optional = true)
    public String adapterMiddleOverride;

    @Argument(fullName = LONG_NAME_NO_5_P_3_P_ADAPTERS, optional = true)
    public boolean no5p3pAdapters;

    @Argument(fullName = LONG_NAME_REVERSE_COMPLEMENT_READ_2, optional = true)
    public boolean reverseComplementRead2 = true;

    @Argument(fullName = LONG_NAME_ILLUMINA, optional = true)
    public boolean illumina;

    @Argument(fullName = LONG_NAME_ILLUMINA_READ_1_LIST, optional = true)
    public List<GATKPath> illuminaRead1List;

    @Argument(fullName = LONG_NAME_ILLUMINA_READ_2_LIST, optional = true)
    public List<GATKPath> illuminaRead2List;

    @Argument(fullName = LONG_NAME_ADAPTER_MIN_ERROR_RATE)
    public double adapterMinErrorRate = 0.2;

    @Argument(fullName = LONG_NAME_ADAPTER_MIN_OVERLAP)
    public int adapterMinOverlap = 10;

    @Argument(fullName = LONG_NAME_QUALITY_CUTOFF)
    public int qualityCutoff = 30;

    // debugging parameters
    @Argument(fullName = FULL_NANE_FASTQ_ASYNC_IO)
    public boolean fastqAsyncIO;

    @Argument(fullName = LONG_NAME_BYPASS_MODE)
    public boolean bypassMode;

    @Argument(fullName = LONG_NAME_NO_OUTPUT)
    public boolean noOutput;

    @Argument(fullName = LONG_NAME_RETURN_FIRST_FOUND_ADAPTER)
    public boolean returnFirstFoundAdapter = true;

    @Argument(fullName = FULL_NAME_LOG_ADAPTERS)
    public LogAdapters logAdapters = LogAdapters.None;

    @Argument(fullName = FULL_NAME_MAX_INPUT_READS)
    public int maxInputReads;

    @Argument(fullName = FULL_NAME_MAX_OUTPUT_READS)
    public int maxOutputReads;

    enum LibraryDirection {
        ThreePrime,
        FivePrime
    }

    enum Chemistry {
        TenX_V2,
        TenX_V3
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
        }

    }
}
