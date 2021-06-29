package org.ultimagen.annotate;

import htsjdk.variant.variantcontext.VariantContext;

public interface VariantContextAnnotator {

    void annotate(final VariantContext vc);
}
