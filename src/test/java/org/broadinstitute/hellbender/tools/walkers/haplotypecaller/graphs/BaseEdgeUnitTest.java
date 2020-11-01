package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class BaseEdgeUnitTest extends GATKBaseTest {
    @DataProvider(name = "EdgeCreationData")
    public Object[][] makeMyDataProvider() {
        List<Object[]> tests = new ArrayList<>();

        // this functionality can be adapted to provide input data for whatever you might want in your data
        for (final int multiplicity : Arrays.asList(1, 2, 3)) {
            for (final boolean isRef : Arrays.asList(true, false)) {
                for (final Boolean forwardStrand : Arrays.asList(null, Boolean.TRUE, Boolean.FALSE)) {
                    tests.add(new Object[]{isRef,
                            forwardStrand == null || forwardStrand ? multiplicity : 0,
                            forwardStrand == null || !forwardStrand ? multiplicity : 0});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "EdgeCreationData")
    public void testBasic(final boolean isRef, final int fMult, final int rMult) {
        final BaseEdge e = new BaseEdge(isRef, fMult, rMult);
        Assert.assertEquals(e.isRef(), isRef);
        Assert.assertEquals(e.getMultiplicity(), fMult + rMult);
        Assert.assertEquals(e.getPruningMultiplicity(), fMult+ rMult);
        Assert.assertEquals(e.getDotLabel(), String.format("%d//%d", fMult, rMult));

        e.toString(); //just check not blowing up

        e.setIsRef(!isRef);
        Assert.assertEquals(e.isRef(), !isRef);

        e.toString();//just check not blowing up

        e.setMultiplicity(fMult + 1, rMult);
        Assert.assertEquals(e.getMultiplicity(), rMult + fMult + 1);
        Assert.assertEquals(e.getPruningMultiplicity(), fMult + 1+ rMult);
        Assert.assertEquals(e.getDotLabel(), String.format("%d//%d", fMult + 1, rMult));

        e.toString();//just check not blowing up

        e.incMultiplicity(1, true);
        Assert.assertEquals(e.getMultiplicity(), rMult + fMult + 2);
        Assert.assertEquals(e.getPruningMultiplicity(), fMult + 1 + rMult + 1);
        Assert.assertEquals(e.getDotLabel(), String.format("%d//%d", fMult + 1, rMult + 1));

        e.toString();//just check not blowing up

        final BaseEdge copy = e.copy();
        Assert.assertEquals(copy.isRef(), e.isRef());
        Assert.assertEquals(copy.getMultiplicity(), e.getMultiplicity());
        Assert.assertEquals(copy.getPruningMultiplicity(), e.getPruningMultiplicity());
        Assert.assertEquals(copy.getDotLabel(), e.getDotLabel());

        e.toString();//just check not blowing up
    }

    @Test(dataProvider = "EdgeCreationData")
    public void testAdd(final boolean isRef, final int fMult,final int rMult) {
        final BaseEdge e1 = new BaseEdge(isRef, fMult,rMult);
        final BaseEdge e2 = new BaseEdge(isRef, fMult,rMult);
        final BaseEdge e3 = e1.add(e2);
        Assert.assertTrue(e1 == e3);//identity
        Assert.assertEquals(e1.isRef(), isRef);
        Assert.assertEquals(e1.getMultiplicity(), 2*(fMult+rMult));
        Assert.assertEquals(e1.getPruningMultiplicity(), 2*(fMult+rMult));
        Assert.assertEquals(e1.getDotLabel(), String.format("%d//%d", 2*fMult, 2*rMult));

        final BaseEdge e4 = new BaseEdge(!isRef, fMult,rMult);
        e1.add(e4);
        Assert.assertTrue(e1.isRef()); //one or the other was ref
    }

    @DataProvider
    Object[][] testAddOrWhichStrand(){
        return new Object[][] {
                {Boolean.TRUE},
                {Boolean.FALSE},
                {null}
        };
    }

    @Test(dataProvider = "testAddOrWhichStrand")
    public void testAddOr(final Boolean forwardStrand) {

        final BaseEdge e1f = new BaseEdge(false, forwardStrand == null || forwardStrand ? 1 : 0, forwardStrand == null || !forwardStrand ? 1 : 0);
        final BaseEdge e2f = new BaseEdge(false, forwardStrand == null || forwardStrand ? 2 : 0, forwardStrand == null || !forwardStrand ? 2 : 0);
        final BaseEdge e1e2 = BaseEdge.makeOREdge(Arrays.asList(e1f, e2f));
        Assert.assertEquals(e1e2.getMultiplicity(), forwardStrand == null ? 6 : 3);
        Assert.assertFalse(e1e2.isRef());

        final BaseEdge e3t = new BaseEdge(true, forwardStrand == null || forwardStrand ? 3 : 0, forwardStrand == null || !forwardStrand ? 3 : 0);
        final BaseEdge e1e2e3 = BaseEdge.makeOREdge(Arrays.asList(e1f, e2f, e3t));
        Assert.assertEquals(e1e2e3.getMultiplicity(), forwardStrand == null ? 12 : 6);
        Assert.assertTrue(e1e2e3.isRef());
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testAddNull() {
        final BaseEdge e = new BaseEdge(false, 1, 2);
        e.add(null);
    }

    @Test(dataProvider = "testAddOrWhichStrand")
    public void testEdgeWeightComparator(final Boolean forwardStrand) {
        final BaseEdge e10 = new BaseEdge(false, forwardStrand == null || forwardStrand ? 10 : 0, forwardStrand == null || !forwardStrand ? 10 : 0);
        final BaseEdge e5 = new BaseEdge(true, forwardStrand == null || forwardStrand ? 5 : 0, forwardStrand == null || !forwardStrand ? 5 : 0);
        final BaseEdge e2 = new BaseEdge(false, forwardStrand == null || forwardStrand ? 2 : 0, forwardStrand == null || !forwardStrand ? 2 : 0);
        final BaseEdge e1 = new BaseEdge(false, forwardStrand == null || forwardStrand ? 1 : 0, forwardStrand == null || !forwardStrand ? 1 : 0);

        final List<BaseEdge> edges = new ArrayList<>(Arrays.asList(e1, e2, e5, e10));
        edges.sort(BaseEdge.EDGE_MULTIPLICITY_ORDER);
        Assert.assertEquals(edges.get(0), e10);
        Assert.assertEquals(edges.get(1), e5);
        Assert.assertEquals(edges.get(2), e2);
        Assert.assertEquals(edges.get(3), e1);
    }
}
