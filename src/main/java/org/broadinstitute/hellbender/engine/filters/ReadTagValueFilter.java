package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.util.Lazy;
import htsjdk.variant.variantcontext.VariantContextUtils;
import org.apache.commons.jexl2.JexlContext;
import org.apache.commons.jexl2.JexlEngine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import org.apache.commons.jexl2.Expression;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
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
            doc = "Look for this tag in read", optional=true)
    public String readFilterTagName = null;

    @Argument(fullName = ReadFilterArgumentDefinitions.READ_FILTER_TAG_COMP,
            doc = "Compare value in tag to this value", optional=true)
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

    @Argument(fullName=ReadFilterArgumentDefinitions.READ_FILTER_EXPRESSION_LONG_NAME, shortName="filter", doc="One or more JEXL expressions used to filter", optional=true)
    public List<String> filterExpressions = new ArrayList<>();


    @Argument(fullName = ReadFilterArgumentDefinitions.READ_FILTER_TAG_OP,
            doc = "Compare value in tag to value with this operator. " +
                    "If T is the value in the tag, OP is the operation provided, " +
                    "and V is the value in read-filter-tag, then the " +
                    "read will pass the filter iff T OP V is true.", optional = true)
    public Operator readFilterTagOp = Operator.EQUAL;

    private Lazy<List<Expression>> jexlExprs = new Lazy<>(() -> {
        List<Expression>        l = new LinkedList<>();
        for ( String expr : filterExpressions )
            l.add(VariantContextUtils.engine.get().createExpression(expr));
        return l;
    });

    private static class GATKReadJexlContext implements JexlContext {

        final private GATKRead        read;

        GATKReadJexlContext(final GATKRead read) {
            this.read = read;
        }

        @Override
        public Object get(String name) {
            return read.getAttributeAsString(name);
        }

        @Override
        public void set(String name, Object value) {
            throw new IllegalArgumentException("setting attributes is not allowed");
        }

        @Override
        public boolean has(String name) {
            return read.hasAttribute(name);
        }
    }

    public ReadTagValueFilter() {
    }

    @Override
    public boolean test(final GATKRead read) {

        // must have either an expression or a name/value arg
        if (  jexlExprs.get().size() == 0 && readFilterTagName == null )
            throw new IllegalArgumentException("must have either an expression or a name/value arg");

        // using jexl?
        if ( jexlExprs.get().size() > 0 ) {

            // loop over expressions. At this point expressions are ANDed
            for ( Expression expr : jexlExprs.get() ) {
                Object v = expr.evaluate(new GATKReadJexlContext(read));
                if (!v.equals(Boolean.TRUE))
                    return false;
            }
        }

        // using name/value operator?
        if ( readFilterTagName != null ) {
            return read.hasAttribute(this.readFilterTagName) &&
                    this.readFilterTagComp != null &&
                    this.readFilterTagOp.comp.apply(read.getAttributeAsFloat(this.readFilterTagName),
                            this.readFilterTagComp);
        } else {
            return true;
        }
    }
}