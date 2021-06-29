package org.ultimagen.annotate;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public abstract class AnnotatorBase implements VariantContextAnnotator {

    // annotation names
    protected final String   A_INDEL_CLASSIFY = "indel_classify";
    protected final String   A_INDEL_LENGTH = "indel_length";
    protected final String   A_HMER_INDEL_LENGTH = "hmer_indel_length";
    protected final String   A_HMER_INDEL_NUC = "hmer_indel_nuc";
    protected final String   A_LEFT_MOTIF = "left_motif";
    protected final String   A_RIGHT_MOTIF = "right_motif";
    protected final String   A_GC_CONTENT = "gc_content";
    protected final String   A_CYCLESKIP_STATUS = "cycleskip_status";

    // additional constants
    protected final String   C_INSERT = "ins";
    protected final String   C_DELETE = "del";
    protected final String   C_CSS_NA = "NA";
    protected final String   C_CSS_CS = "cycle-skip";
    protected final String   C_CSS_PCS = "possible-cycle-skip";
    protected final String   C_CSS_NS = "non-skip";


    // bean fields
    private ReferenceSequenceFile referenceSequenceFile;

    public AnnotatorBase(final ReferenceSequenceFile referenceSequenceFile) {
        this.referenceSequenceFile = referenceSequenceFile;
    }

    // get nucleoid from reference
    protected byte getReferenceNucleoid(final String contig, final int start) {
        return referenceSequenceFile.getSubsequenceAt(contig, start, start).getBases()[0];
    }

    // get motif from reference
    protected String getReferenceMotif(final String contig, final int start, final int end) {
        return referenceSequenceFile.getSubsequenceAt(contig, start, end).getBaseString();
    }


    // add a single 'cooked' annotation
    protected void annotate(final VariantContext vc, final String name, final Object value) {
        if ( value != null ) {
            vc.getAttributes().put(name, value);
        }
    }
}
