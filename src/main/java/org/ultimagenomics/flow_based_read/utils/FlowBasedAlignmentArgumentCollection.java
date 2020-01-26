package org.ultimagenomics.flow_based_read.utils;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;

import java.io.Serializable;

public class FlowBasedAlignmentArgumentCollection implements Serializable {
    private static final long serialVersionUID = -3405667551696313040L;

    private static final String PROBABILITY_RATIO_THRESHOLD_LONG_NAME = "flow-probability-threshold";
    private static final String REMOVE_LONGER_THAN_ONE_INDELS = "flow-remove-non-single-base-pair-indels";
    private static final String REMOVE_ONE_TO_ZERO_PROBS = "flow-remove-one-zero-probs";
    private static final String NUMBER_OF_POSSIBLE_PROBS = "flow-quantization-bins";
    private static final String FILLING_VALUE = "flow-fill-empty-bins-value";
    private static final String SYMMETRIC_INDELS = "flow-symmetric-indel-probs";
    private static final String REPORT_INS_OR_DEL = "flow-report-insertion-or-deletion";
    private static final String DISALLOW_LARGER_PROBS = "flow-disallow-probs-larger-than-call";
    private static final String LUMP_PROBS = "flow-lump-probs";
    private static final String PROB_SF = "flow-probability-scaling-factor";
    private static final String RETAIN_MAX_N_PROBS_BASE = "flow-retain-max-n-probs-base-format";

    private static final double DEFAULT_RATIO_THRESHOLD = 0.003;
    private static final double DEFAULT_FILLING_VALUE = 0.001;
    private static final boolean DEFAULT_REMOVE_LONGER_INDELS = false;
    private static final boolean DEFAULT_REMOVE_ONE_TO_ZERO = false;
    private static final boolean DEFAULT_SYMMETRIC_INDELS = false;
    private static final int DEFAULT_QUANTIZATION = 121;
    private static final boolean DEFAULT_ONLY_INS_OR_DEL = false;
    private static final boolean DEFAULT_DISALLOW_LARGER_PROBS = false;
    private static final boolean DEFAULT_LUMP_PROBS = false;
    private static final boolean DEFAULT_RETAIN_MAX_N_PROBS = false;
    private static final int DEFAULT_PROB_SCALING_FACTOR = 10;

    @Advanced
    @Argument(fullName = PROBABILITY_RATIO_THRESHOLD_LONG_NAME, doc = "Lowest probability ratio to be used as an option", optional = true)
    public double probability_ratio_threshold = DEFAULT_RATIO_THRESHOLD;

    @Advanced
    @Argument(fullName = REMOVE_LONGER_THAN_ONE_INDELS, doc = "Should the probabilities of more then 1 indel be used", optional = true)
    public boolean remove_longer_than_one_indels = DEFAULT_REMOVE_LONGER_INDELS;

    @Advanced
    @Argument(fullName = REMOVE_ONE_TO_ZERO_PROBS, doc = "Remove probabilities of basecall of zero from non-zero genome", optional = true)
    public boolean remove_one_to_zero_probs = DEFAULT_REMOVE_ONE_TO_ZERO;

    @Advanced
    @Argument(fullName = NUMBER_OF_POSSIBLE_PROBS, doc = "Probability quantization", optional = true)
    public int probability_quantization = DEFAULT_QUANTIZATION;

    @Advanced
    @Argument(fullName = FILLING_VALUE, doc = "Value to fill the zeros of the matrix with", optional=true)
    public double filling_value = DEFAULT_FILLING_VALUE;

    @Advanced
    @Argument(fullName = SYMMETRIC_INDELS, doc = "Should indel probabilities be symmetric in flow", optional=true)
    public boolean symmetric_indels = DEFAULT_SYMMETRIC_INDELS;

    @Advanced
    @Argument(fullName = REPORT_INS_OR_DEL, doc = "Report either insertion or deletion, probability, not both", optional=true)
    public boolean only_ins_or_del = DEFAULT_ONLY_INS_OR_DEL;

    @Advanced
    @Argument(fullName = DISALLOW_LARGER_PROBS, doc = "Cap probabilities of error to 1 relative to base call", optional=true)
    public boolean disallow_larger_probs = DEFAULT_DISALLOW_LARGER_PROBS;

    @Advanced
    @Argument(fullName = LUMP_PROBS, doc = "Should all probabilities of insertion or deletion in the flow be combined together", optional=true)
    public boolean lump_probs = DEFAULT_LUMP_PROBS;

    @Advanced
    @Argument(fullName = RETAIN_MAX_N_PROBS_BASE, doc = "Keep only hmer/2 probabilities (like in base format)", optional=true)
    public boolean retainMaxNProbs = DEFAULT_RETAIN_MAX_N_PROBS;

    @Advanced
    @Argument(fullName = PROB_SF, doc = "probability scaling factor for (phred=10)", optional=true)
    public int probability_scaling_factor = DEFAULT_PROB_SCALING_FACTOR;

    public FlowBasedAlignmentArgumentCollection() {}

}
