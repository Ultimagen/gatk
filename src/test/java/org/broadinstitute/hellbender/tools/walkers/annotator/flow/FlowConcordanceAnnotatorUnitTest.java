package org.broadinstitute.hellbender.tools.walkers.annotator.flow;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotator;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

public class FlowConcordanceAnnotatorUnitTest {

    @Test
    public void testBasic() {

        final String[][]        testData = {
                // order:
                // refbases, altAllele (include a space before and after refAllele
                // indel-class, indel-length, hmer-indel-lenfth, hmer-indel-nuc
                // left-motif, right-motif, gc-content, cycleskip-status
                // value that starts with "!" means ignore
                {
                        // a simple SNP
                        "GTATC A ACATCGGA", "C",
                        "NA", null, null, null, "GTATC", "ACATC", "0.3", "non-skip"
                },
                {
                        // a possible-cycle-skip SNP
                        "GTATC A TCATCGGA", "C",
                        "NA", null, null, null, "GTATC", "TCATC", "0.3", "possible-cycle-skip"
                },
                {
                        // a cycle-skip SNP
                        "GTATC A ACATCGGA", "T",
                        "NA", null, null, null, "GTATC", "ACATC", "0.3", "cycle-skip"
                },
                {
                        // not hmer indel
                        "TATCT CA TTGACCAA", "C",
                        "del", "1", null, null, "ATCTC", "TTGAC", "0.3", "NA"
                },
                {
                        // del hmer indel
                        "TATCTC AT TGACCAA", "A",
                        "del", "1", "2", "T", "TCTCA", "GACCA", "0.4", "NA"
                },
                {
                        // ins hmer indel
                        "TATCT C ATTGACCAA", "CA",
                        "ins", "1", "1", "A", "ATCTC", "TTGAC", "0.3", "NA"
                }
        };

        // setup flow order
        VariantAnnotator.flowOrder = "TGCA";

        // should be in same order as test data!!!!
        final List<String>      expectedAttrs = (new FlowConcordanceAnnotator()).getKeyNames();

        // loop on test data
        for ( String[] data : testData ) {

            // prepare
            final int        refAlleleStart = data[0].indexOf(' ');
            final int        refAlleleEnd = data[0].indexOf(' ', refAlleleStart + 1);
            final String     refAllele = data[0].substring(refAlleleStart + 1, refAlleleEnd);
            final ReferenceContext ref = buildReferenceContext(data[0].replace(" ", ""), refAlleleStart + 1, refAlleleEnd - 1);
            final VariantContext vc = buildVariantContext(ref, refAllele, data[1]);
            String          msg = "on " + StringUtils.join(data, " ");

            // invoke
            final Map<String, Object> attrs = FlowConcordanceAnnotator.annotateForTesting(ref, vc);
            Assert.assertNotNull(attrs, msg);

            // check that all expected attributes are there
            for ( int n = 0 ; n < 8 ; n++ ) {
                String       key = expectedAttrs.get(n);
                String       elem = data[2+n];
                String       keyMsg = "on " + key + " " + msg;
                if ( elem != null && elem.charAt(0) != '!' ) {
                    Object v = attrs.get(key);
                    if (v instanceof List) {
                        v = StringUtils.join((List) v, ",");
                    }
                    Assert.assertEquals(v.toString(), elem, keyMsg);
                } else if ( elem == null ) {
                    Assert.assertFalse(attrs.containsKey(key), keyMsg);
                }
            }
        }
    }

    private ReferenceContext buildReferenceContext(String refBases, int start, int stop) {

        // note that locations here are 1 based
        final String                insLoc = "chr1";
        final SimpleInterval        interval = new SimpleInterval(insLoc, start, stop);
        final byte[]                refBytes = refBases.getBytes();
        final SimpleInterval        interval1 = new SimpleInterval(insLoc, 1, refBytes.length);
        final ReferenceBases        ref1 = new ReferenceBases(refBytes, interval1);
        final SAMSequenceDictionary dict = new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord(insLoc, refBytes.length)));
        final ReferenceContext      ref = new ReferenceContext(ReferenceDataSource.of(ref1, dict), interval, start - 1, 20);

        return ref;
    }

    private VariantContext buildVariantContext(ReferenceContext ref, String refBases, String altBases) {

        final Allele refAllele = Allele.create(refBases, true);
        final Allele altAllele = Allele.create(altBases, false);
        final VariantContext vc = new VariantContextBuilder("foo",
                ref.getContig(), ref.getStart(), ref.getEnd(),
                Arrays.asList(refAllele, altAllele)).make();

        return vc;
    }
}