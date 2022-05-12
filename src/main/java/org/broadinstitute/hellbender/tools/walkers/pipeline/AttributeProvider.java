package org.broadinstitute.hellbender.tools.walkers.pipeline;

public interface AttributeProvider {

    boolean hasAttribute(final String attributeName);
    int getAttributeAsInteger(final String attributeName);
}
