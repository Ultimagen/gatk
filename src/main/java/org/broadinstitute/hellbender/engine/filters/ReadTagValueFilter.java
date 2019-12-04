package org.broadinstitute.hellbender.engine.filters;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.function.BiFunction;

/**
 * Keep only reads that contain a tag with a value that agrees with parameters as specified.
 *
 */
@DocumentedFeature(groupName = HelpConstants.DOC_CAT_READFILTERS, groupSummary = HelpConstants.DOC_CAT_READFILTERS_SUMMARY,
        summary = "Keep only reads that contains a tag with value that agrees with parameters")
public final class ReadTagValueFilter extends ReadFilter {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = ReadFilterArgumentDefinitions.READ_FILTER_TAG,
            doc = "Look for this tag in read")
    public String readFilterTagName = null;

    @Argument(fullName = ReadFilterArgumentDefinitions.READ_FILTER_TAG_COMP,
            doc = "Compare value in tag to this value")
    public Float readFilterTagComp = 0F;

    public enum Operator {
        LESS((Float x, Float y) -> x < y),
        LESS_OR_EQUAL((Float x, Float y) -> x <= y),
        GREATER((Float x, Float y) -> x > y),
        GREATER_OR_EQUAL((Float x, Float y) -> x >= y),
        EQUAL(Float::equals),
        NOT_EQUAL((Float x, Float y) -> !x.equals(y));

        BiFunction<Float, Float, Boolean> comp;

        Operator(BiFunction<Float, Float, Boolean> comp) {
            this.comp = comp;
        }
    }

    @Argument(fullName = ReadFilterArgumentDefinitions.READ_FILTER_TAG_OP,
            doc = "Compare value in tag to value with this operator. " +
                    "If T is the value in the tag, OP is the operation provided, " +
                    "and V is the value in read-filter-tag, then the " +
                    "read will pass the filter iff T OP V is true.")
    public Operator readFilterTagOp = Operator.EQUAL;

    public ReadTagValueFilter() {}

    @Override
    public boolean test(final GATKRead read) {
        return this.readFilterTagName != null &&
                read.hasAttribute(this.readFilterTagName) &&
                this.readFilterTagComp != null &&
                this.readFilterTagOp.comp.apply(read.getAttributeAsFloat(this.readFilterTagName),
                        this.readFilterTagComp);
    }
}
