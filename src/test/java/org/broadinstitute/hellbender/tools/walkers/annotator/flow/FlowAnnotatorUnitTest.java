package org.broadinstitute.hellbender.tools.walkers.annotator.flow;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.util.*;

public class FlowAnnotatorUnitTest {

    private FlowAnnotatorBase[]     allAnnotators = {
            new IndelClassify(),
            new IndelLength(),
            new HmerIndelLength(),
            new HmerIndelNuc(),
            new LeftMotif(),
            new RightMotif(),
            new GcContent(),
            new CycleSkipStatus()
    };

    final static String         FLOW_ORDER = "TGCA";
    final static int            RANDOM_TEST_DATA_ENTRY_COUNT = 1000;
    final static int            RANDOM_TEST_DATA_MIN_REF_LENGTH = 5;
    final static int            RANDOM_TEST_DATA_MAX_REF_LENGTH = 10;
    final static int            RANDOM_TEST_DATA_MIN_ALLELE_LENGTH = 1;
    final static int            RANDOM_TEST_DATA_MAX_ALLELE_LENGTH = 2;

    @DataProvider(name = "testData")
    public Object[][] getTestData() {

        final Object[][]        testData = {
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

        return testData;
    }

    @DataProvider(name = "randomTestData")
    public Object[][] getRandomTestData() {

        // this test data is designed to increase confidence in the code no crashing
        final Object[][]        testData = new Object[RANDOM_TEST_DATA_ENTRY_COUNT][];
        final Random            rand = new Random(0);
        for ( int i = 0 ; i < testData.length ; i++ ) {

            // random entries consist of:
            // [0] reference
            // [1] allele
            // [2] index of annotation to call
            testData[i] = new Object[3];
            final String    alleleOnRef = generateRandomSequence(rand, FLOW_ORDER, RANDOM_TEST_DATA_MIN_ALLELE_LENGTH, RANDOM_TEST_DATA_MAX_ALLELE_LENGTH, null);
            final String    allele = generateRandomSequence(rand, FLOW_ORDER, RANDOM_TEST_DATA_MIN_ALLELE_LENGTH, RANDOM_TEST_DATA_MAX_ALLELE_LENGTH, alleleOnRef);
            testData[i][0] =
                    generateRandomSequence(rand, FLOW_ORDER, RANDOM_TEST_DATA_MIN_REF_LENGTH, RANDOM_TEST_DATA_MAX_REF_LENGTH, null)
                    + " "
                    + alleleOnRef
                    + " "
                    + generateRandomSequence(rand, FLOW_ORDER, RANDOM_TEST_DATA_MIN_REF_LENGTH, RANDOM_TEST_DATA_MAX_REF_LENGTH, null);
            testData[i][1] = allele;
            testData[i][2] = rand.nextInt(allAnnotators.length);
        }

        return testData;
    }

    private String generateRandomSequence(final Random rand, final String flowOrder,
                                          final int minLength, final int maxLength,
                                          final String exceptFor) {

        final int           length = minLength + rand.nextInt(maxLength - minLength);
        do {
            final StringBuilder sb = new StringBuilder();
            for (int i = 0; i < length; i++)
                sb.append(flowOrder.charAt(rand.nextInt(flowOrder.length())));

            if ( exceptFor == null || !sb.toString().equals(exceptFor) )
                return sb.toString();
        } while (true);
    }

    @Test(dataProvider = "testData")
    public void testBasic(Object[] testData) {

        // should be in same order as test data!!!!
        final List<String>      expectedAttrs = allKeys();

        // prepare as array of Strings
        String[]        data = new String[testData.length];
        for ( int n = 0 ; n < data.length ; n++ ) {
            if (testData[n] != null) {
                data[n] = testData[n].toString();
            }
        }

        // prepare
        final int        refAlleleStart = data[0].indexOf(' ');
        final int        refAlleleEnd = data[0].indexOf(' ', refAlleleStart + 1);
        final String     refAllele = data[0].substring(refAlleleStart + 1, refAlleleEnd);
        final ReferenceContext ref = buildReferenceContext(data[0].replace(" ", ""), refAlleleStart + 1, refAlleleEnd - 1);
        final VariantContext vc = buildVariantContext(ref, refAllele, data[1]);
        String          msg = "on " + StringUtils.join(data, " ");

        // invoke
        final Map<String, Object> attrs = allAnnotate(ref, vc);
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

    @Test(dataProvider = "randomTestData")
    public void testRandom(Object[] testData) {

        // should be in same order as test data!!!!
        final List<String>      expectedAttrs = allKeys();

        // prepare specific types
        final String    data0 = testData[0].toString();
        final String    data1 = testData[1].toString();
        final int       data2 = Integer.parseInt(testData[2].toString());

        // prepare
        final int        refAlleleStart = data0.indexOf(' ');
        final int        refAlleleEnd = data0.indexOf(' ', refAlleleStart + 1);
        final String     refAllele = data0.substring(refAlleleStart + 1, refAlleleEnd);
        final ReferenceContext ref = buildReferenceContext(data0.replace(" ", ""), refAlleleStart + 1, refAlleleEnd - 1);
        final VariantContext vc = buildVariantContext(ref, refAllele, data1);
        String          msg = "on " + StringUtils.join(testData, " ");

        // invoke
        final Map<String, Object> attrs = allAnnotators[data2].annotate(ref, vc, null);
        Assert.assertNotNull(attrs, msg);
    }

    private Map<String, Object> allAnnotate(final ReferenceContext ref, final VariantContext vc) {

        final Map<String, Object>     attrs = new LinkedHashMap<>();

        for ( FlowAnnotatorBase a : allAnnotators ) {
            attrs.putAll(a.annotate(ref, vc, null));
        }

        return attrs;
    }

    private List<String> allKeys() {

        List<String>     keys = new LinkedList<>();

        for ( FlowAnnotatorBase a : allAnnotators ) {
            keys.addAll(a.getKeyNames());
            a.setFlowOrder(Collections.singletonList(FLOW_ORDER));
        }

        return keys;

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