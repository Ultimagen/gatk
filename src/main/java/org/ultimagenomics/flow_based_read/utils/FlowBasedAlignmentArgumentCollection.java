package org.ultimagenomics.flow_based_read.utils;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

public class FlowBasedAlignmentArgumentCollection implements Serializable {
    private static final long serialVersionUID = -3405667551696313040L;

    public static final String PROBABILITY_RATIO_THRESHOLD_LONG_NAME = "flow-probability-threshold";
    public static final String REMOVE_LONGER_THAN_ONE_INDELS = "flow-remove-non-single-base-pair-indels";
    public static final String REMOVE_ONE_TO_ZERO_PROBS = "flow-remove-one-zero-probs";
    public static final String NUMBER_OF_POSSIBLE_PROBS = "flow-quantization-bins";
    public static final String FILLING_VALUE="flow-fill-empty-bins-value";

    public static final double DEFAULT_RATIO_THRESHOLD = 0.003;
    public static final double DEFAULT_FILLING_VALUE = 0.001;
    public static final boolean DEFAULT_REMOVE_LONGER_INDELS = false;
    public static final boolean DEFAULT_REMOVE_ONE_TO_ZERO = false;
    public static final int DEFAULT_QUANTIZATION = 121;

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

    public FlowBasedAlignmentArgumentCollection() {}
}




