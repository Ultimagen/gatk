package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.*;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

public class ReadTagValueFilterUnitTest extends GATKBaseTest {

    private static final int CHR_COUNT = 1;
    private static final int CHR_START = 1;
    private static final int CHR_SIZE = 1000;
    private static final int GROUP_COUNT = 5;

    private final SAMFileHeader header= ArtificialReadUtils.createArtificialSamHeaderWithGroups(CHR_COUNT, CHR_START, CHR_SIZE, GROUP_COUNT);

    private GATKRead buildSAMRead(final String cigarString, Object[] attrsNameAndValue) {
        final Cigar cigar = TextCigarCodec.decode(cigarString);
        final SAMRecord samRecord = ArtificialReadUtils.createArtificialSAMRecord(header, cigar);
        for ( int i = 0 ; i < attrsNameAndValue.length ; i += 2 )
            samRecord.setAttribute(attrsNameAndValue[i].toString(), attrsNameAndValue[i+1]);
        return new SAMRecordToGATKReadAdapter(samRecord);
    }

    @Test(dataProvider= "ReadTagValueFilterDataProvider")
    public void testReadTagValueFilter(final String cigarString,
                                      final Object[] attrsNameAndValue,
                                      final String tagName,
                                      final Float tagValue,
                                      final ReadTagValueFilter.Operator tagOp,
                                      final String[] jexlExpr,
                                      final boolean expectedResult) {

        final ReadTagValueFilter filter;

        // test different constructors here as well
        if ( tagName != null && jexlExpr.length == 0 ) {
            filter = new ReadTagValueFilter(tagName, tagValue, tagOp);
        } else if ( tagName == null && jexlExpr.length == 1 ) {
            filter = new ReadTagValueFilter(jexlExpr[0]);
        } else if ( tagName == null && jexlExpr.length > 1 ) {
            filter = new ReadTagValueFilter(Arrays.asList(jexlExpr));
        } else {
            filter = new ReadTagValueFilter();
            filter.readFilterTagName = tagName;
            filter.readFilterTagComp = tagValue;
            filter.readFilterTagOp = tagOp;
            filter.filterExpressions = Arrays.asList(jexlExpr);
        }

        final GATKRead read = buildSAMRead(cigarString, attrsNameAndValue);
        Assert.assertEquals(filter.test(read), expectedResult, cigarString);
    }

    @DataProvider(name = "ReadTagValueFilterDataProvider")
    public Iterator<Object[]> readTagValueFilterDataProvider() {
        final List<Object[]> result = new LinkedList<>();

        result.add(new Object[] {
                "100M",                             // vigat
                new Object[] {"TM", 1.0f},           // attributes
                "TM",                               // tagname
                1.0f,                               // tagvalue
                ReadTagValueFilter.Operator.EQUAL,  // tagop
                new String[] {},                    // jexl expressions
                Boolean.TRUE                        // expected
        });

        result.add(new Object[] {
                "100M",                             // vigat
                new Object[] {"TM", 1.0f},           // attributes
                "TM",                               // tagname
                1.0f,                               // tagvalue
                ReadTagValueFilter.Operator.LESS,  // tagop
                new String[] {},                    // jexl expressions
                Boolean.FALSE                       // expected
        });

        result.add(new Object[] {
                "100M",                             // vigat
                new Object[] {"TM", 1.0f},           // attributes
                null,                               // tagname
                null,                               // tagvalue
                null,                               // tagop
                new String[] {"TM == 1.0"},         // jexl expressions
                Boolean.TRUE                       // expected
        });

        result.add(new Object[] {
                "100M",                             // vigat
                new Object[] {"TM", 1.0f},           // attributes
                null,                               // tagname
                null,                               // tagvalue
                null,                               // tagop
                new String[] {"TM < 1.0"},         // jexl expressions
                Boolean.FALSE                       // expected
        });

        result.add(new Object[] {
                "100M",                             // vigat
                new Object[] {"TM", 1.0f},           // attributes
                "TM",                               // tagname
                1.0f,                               // tagvalue
                ReadTagValueFilter.Operator.LESS,  // tagop
                new String[] {"TM == 1.0"},         // jexl expressions
                Boolean.FALSE                       // expected
        });
        return result.iterator();
    }


}
