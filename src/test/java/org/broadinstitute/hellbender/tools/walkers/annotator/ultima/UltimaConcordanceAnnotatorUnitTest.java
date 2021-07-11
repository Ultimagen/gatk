package org.broadinstitute.hellbender.tools.walkers.annotator.ultima;

import com.google.errorprone.annotations.Var;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.tools.walkers.annotator.CountNs;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

public class UltimaConcordanceAnnotatorUnitTest {

    @Test
    public void testBasic() {

        final Object[][]        testData = {
                {
                    // for now only a basic test isd defined. more could be added below using the same format
                    "GTATCATCATCGGA", 6, 6, "A", "AATC",                // refbases, start, stop, refAllele, altAllele
                    UltimaConcordanceAnnotator.C_INSERT, 3, 1, "A",     // indel-classify, indel-lenfth, hmer-indel-length, hmer-indel-nuc
                    "GTATC", "TCATC", (float)0.3, UltimaConcordanceAnnotator.C_CSS_NA // left-motif, right-motif, gc-content, cycleskip-status
                }
        };

        // should be in same order as test data!!!!
        final List<String>      expectedAttrs = (new UltimaConcordanceAnnotator()).getKeyNames();

        // loop on test data
        for ( Object[] data : testData ) {

            // prepare
            final ReferenceContext ref = buildReferenceContext(data[0].toString(), Integer.parseInt(data[1].toString()), Integer.parseInt(data[2].toString()));
            final VariantContext vc = buildVariantContext(ref, data[3].toString(), data[4].toString());

            // invoke
            final Map<String, Object> attrs = UltimaConcordanceAnnotator.annotateForTesting(ref, vc);
            Assert.assertNotNull(attrs);

            // check that all expected attributes are there
            for ( int n = 0 ; n < 8 ; n++ ) {
                Assert.assertEquals(attrs.get(expectedAttrs.get(n)), data[5+n]);
            }
        }
    }

    private ReferenceContext buildReferenceContext(String refBases, int start, int stop) {

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