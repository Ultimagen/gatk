package org.broadinstitute.hellbender.tools.walkers.annotator.ultima;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.GenotypeAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardMutectAnnotation;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.util.List;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

public abstract class UltimaAnnotatorBase extends GenotypeAnnotation implements StandardMutectAnnotation {

    // additional constants
    protected final String   C_INSERT = "ins";
    protected final String   C_DELETE = "del";
    protected final String   C_CSS_NA = "NA";
    protected final String   C_CSS_CS = "cycle-skip";
    protected final String   C_CSS_PCS = "possible-cycle-skip";
    protected final String   C_CSS_NS = "non-skip";

    @Override
    public void annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        Utils.nonNull(vc);
        Utils.nonNull(g);
        Utils.nonNull(gb);

        annotate(vc, (name, value) -> {
            if ( value != null )
                gb.attribute(name, value);
        });
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return getKeyNames().stream().map(s -> GATKVCFHeaderLines.getFormatLine(s)).collect(Collectors.toList());
    }

    public abstract void annotate(final VariantContext vc, BiConsumer<String, Object> annotator);
}
