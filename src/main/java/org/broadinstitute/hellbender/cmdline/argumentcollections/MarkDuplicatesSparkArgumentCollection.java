package org.broadinstitute.hellbender.cmdline.argumentcollections;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.read.markduplicates.MarkDuplicatesScoringStrategy;
import picard.sam.markduplicates.MarkDuplicates;

import java.io.Serializable;
import java.util.List;


/**
 * An argument collection for use with tools that mark optical
 * duplicates.
 */
public final class MarkDuplicatesSparkArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final String DO_NOT_MARK_UNMAPPED_MATES_LONG_NAME = "do-not-mark-unmapped-mates";
    public static final String DUPLICATE_TAGGING_POLICY_LONG_NAME = "duplicate-tagging-policy";
    public static final String REMOVE_ALL_DUPLICATE_READS = "remove-all-duplicates";
    public static final String REMOVE_SEQUENCING_DUPLICATE_READS = "remove-sequencing-duplicates";

    @Argument(shortName = StandardArgumentDefinitions.DUPLICATE_SCORING_STRATEGY_SHORT_NAME, fullName = StandardArgumentDefinitions.DUPLICATE_SCORING_STRATEGY_LONG_NAME, doc = "The scoring strategy for choosing the non-duplicate among candidates.")
    public MarkDuplicatesScoringStrategy duplicatesScoringStrategy = MarkDuplicatesScoringStrategy.SUM_OF_BASE_QUALITIES;

    @Argument(fullName = MarkDuplicatesSparkArgumentCollection.DO_NOT_MARK_UNMAPPED_MATES_LONG_NAME, doc = "Enabling this option will mean unmapped mates of duplicate marked reads will not be marked as duplicates.")
    public boolean dontMarkUnmappedMates = false;

    // added params
    @Argument(doc = "tagging policy settings. values: DontTag, OpticalOnly, All. default: DontTag")
    public MarkDuplicates.DuplicateTaggingPolicy taggingPolicy = MarkDuplicates.DuplicateTaggingPolicy.DontTag;

    @Argument(doc = "removeSequencingDuplicates. defualt: false")
    public boolean removeSequencingDuplicates = false;

    @Argument(doc = "removeAllDuplicates. default: false")
    public boolean removeAllDuplicates = false;

    @Argument(doc = "Use specific quality summing strategy for flow based reads. The strategy ensures that the same " +
            "(and correct) quality value is used for all bases of the same homopolymer. Default false.")
    public boolean FLOW_QUALITY_SUM_STRATEGY = false;

    @Argument(doc = "Make end location of read be significant when considering duplicates, " +
            "in addition to the start location, which is always significant. Default false.")
    public boolean FLOW_END_LOCATION_SIGNIFICANT = false;

    @Argument(doc = "Maximal number of bases of reads ends difference that is marked as match. Default 0.")
    public int ENDS_READ_UNCERTAINTY = 0;

    @Argument(doc = "Use clipped, rather than unclipped, when considering duplicates. Default false.")
    public boolean FLOW_USE_CLIPPED_LOCATIONS = false;

    @Argument(doc = "Skip first N flows, when considering duplicates. Default 0.")
    public int FLOW_SKIP_START_HOMOPOLYMERS = 0;

    @Argument(doc = "Emit additional debugging info specific to ultima flow. Default false.")
    public boolean DEBUG_ULTIMA_DUPS = false;

    @Argument(doc = "Emit additional debugging info specific to ultima flow: read name. Default null", optional = true)
    public List<String> DEBUG_ULTIMA_READ_NAME = null;

}