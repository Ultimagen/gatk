package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.base.Strings;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingEngine;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrReferenceAnalyzer;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanJavaAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Unit tests for {@link HaplotypeCallerGenotypingEngine}.
 */
public final class HaplotypeCallerGenotypingEngineUnitTest extends GATKBaseTest {
    static String DRAGEN_GATK_BQDFRD_TEST_BAM_SRA = "/Users/emeryj/hellbender/DRAGENMatlab/frdbqd/SRA056922_hs37d5_xmapq.bam";

    private class BasicGenotypingTestProvider extends TestDataProvider {
        byte[] ref;
        byte[] hap;
        Map<Integer,Byte> expected;

        public BasicGenotypingTestProvider(String refString, String hapString, Map<Integer, Byte> expected) {
            super(BasicGenotypingTestProvider.class, String.format("Haplotype to VCF test: ref = %s, alignment = %s", refString,hapString));
            ref = refString.getBytes();
            hap = hapString.getBytes();
            this.expected = expected;
        }
        
        public Map<Integer, Event> calcAlignment() {
            final SmithWatermanAlignment alignment = SmithWatermanJavaAligner.getInstance().align(ref, hap, new SWParameters(3, -1, -4, -1), SWOverhangStrategy.SOFTCLIP);
            final Haplotype h = new Haplotype(hap, false, alignment.getAlignmentOffset(), alignment.getCigar());
            return EventMap.fromHaplotype(h, ref, new SimpleInterval("4", 1, 1 + ref.length), 1);
        }

        public String toString() {
            return "REF:" + new String(ref) + ",ALT:" + new String(hap);
        }
    }


    @DataProvider(name = "BasicGenotypingTestProvider")
    public Object[][] makeBasicGenotypingTests() {

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<>();
            map.put(1 + contextSize, (byte)'M');
            final String context = Strings.repeat("G", contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGACTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<>();
            map.put(2 + contextSize, (byte)'M');
            map.put(21 + contextSize, (byte)'M');
            final String context = Strings.repeat("G", contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG", "ATCTCGCATCGCGAGCATCGCCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'I');
            final String context = Strings.repeat("G", contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGACACTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'D');
            final String context = Strings.repeat("G", contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<>();
            map.put(1, (byte)'M');
            map.put(20, (byte)'D');
            final String context = Strings.repeat("G", contextSize);
            new BasicGenotypingTestProvider("AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCGATAG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<>();
            map.put(2 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'I');
            map.put(30 + contextSize, (byte)'D');
            final String context = Strings.repeat("G", contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "ACCTCGCATCGCGAGCATCGTTACTAGCCGATG", map);
        }

        for( int contextSize : new int[]{0,1,5,9,24,36} ) {
            Map<Integer, Byte> map = new HashMap<>();
            map.put(1 + contextSize, (byte)'M');
            map.put(20 + contextSize, (byte)'D');
            map.put(28 + contextSize, (byte)'M');
            final String context = Strings.repeat("G", contextSize);
            new BasicGenotypingTestProvider(context + "AGCTCGCATCGCGAGCATCGACTAGCCGATAG" + context, "CGCTCGCATCGCGAGCATCGCTAGCCCATAG", map);
        }

        return BasicGenotypingTestProvider.getTests(BasicGenotypingTestProvider.class);
    }
    
    @Test(dataProvider = "BasicGenotypingTestProvider", enabled = true)
    public void testHaplotypeToVCF(BasicGenotypingTestProvider cfg) {
        Map<Integer, Event> calculatedMap = cfg.calcAlignment();
        Map<Integer,Byte> expectedMap = cfg.expected;
        logger.warn(String.format("Test: %s", cfg.toString()));
        if(!compareVCMaps(calculatedMap, expectedMap)) {
            logger.warn("calc map = " + calculatedMap);
            logger.warn("expected map = " + expectedMap);
        }
        Assert.assertTrue(compareVCMaps(calculatedMap, expectedMap),"" + cfg);
    }

    @Test(dataProvider="AddMiscellaneousDataProvider", enabled=false)
    public void testAddMiscellaneousAllele(final String readBases, final int readOffset,
                                           final String ref, final int refOffset,
                                           final String referenceAllele, final String[] alternatives, final double[] likelihoods, final double[] expected) {
        final byte baseQual = (byte)30;

        final byte[] baseQuals = Utils.dupBytes(baseQual, readBases.length());
        final GATKRead read = ArtificialReadUtils.createArtificialRead(readBases.getBytes(), baseQuals, readBases.length() + "M");
        final Locatable loc = new SimpleInterval("20",refOffset,refOffset);
        final ReadPileup pileup = new ReadPileup(loc, Collections.singletonList(read),readOffset);
        final VariantContextBuilder vcb = new VariantContextBuilder();
        final GenotypeBuilder gb = new GenotypeBuilder();
        final List<String> alleleStrings = new ArrayList<>( 1 + alternatives.length);
        alleleStrings.add(referenceAllele);
        alleleStrings.addAll(Arrays.asList(alternatives));

        gb.AD(new int[] { 1 });
        gb.DP(1);
        gb.PL(likelihoods);

        vcb.alleles(alleleStrings);
        vcb.loc("20",refOffset,refOffset + referenceAllele.length() -1);

        vcb.genotypes(gb.make());

        final VariantContext vc = vcb.make();

        final VariantContext updatedVc = null; // GenotypingEngine.addMiscellaneousAllele(vc,pileup,ref.getBytes(),0);
        final GenotypeLikelihoods updatedLikelihoods = updatedVc.getGenotype(0).getLikelihoods();
        Assert.assertEquals(updatedLikelihoods.getAsVector().length, expected.length);
        final double[] updatedLikelihoodsArray = updatedVc.getGenotype(0).getLikelihoods().getAsVector();
        for (int i = 0; i < updatedLikelihoodsArray.length; i++) {
            Assert.assertEquals(updatedLikelihoodsArray[i],expected[i],0.0001);
        }
        Allele altAllele = null;
        for (final Allele allele : updatedVc.getAlleles())
            if (allele.isSymbolic() && allele.getBaseString().equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE_NAME))
                altAllele = allele;
        Assert.assertNotNull(altAllele);
    }

    @DataProvider(name="AddMiscellaneousDataProvider")
    public Iterator<Object[]> addMiscellaneousAlleleDataProvider() {
        return Arrays.asList(ADD_MISCELLANEOUS_ALLELE_DATA).iterator();
    }

    private static final double MATCH_LnLK = QualityUtils.qualToProbLog10((byte)30);
    private static final double MISS_LnLK = QualityUtils.qualToErrorProbLog10((byte)30);

    private static final Object[][] ADD_MISCELLANEOUS_ALLELE_DATA = new Object[][] {
            new Object[] {"ACTG", 0,"ACTGTGAGTATTCC",0,"A",new String[]{}, new double[] {MATCH_LnLK * MATCH_LnLK}, 6,
                    new double[] {MATCH_LnLK * MATCH_LnLK,MATCH_LnLK * MISS_LnLK, MISS_LnLK * MISS_LnLK}}
    };

    /**
     * Private function to compare Map of VCs, it only checks the types and start locations of the VariantContext
     */
    private boolean compareVCMaps(Map<Integer, Event> calc, Map<Integer, Byte> expected) {
        if( !calc.keySet().equals(expected.keySet()) ) { return false; } // sanity check
        for( Integer loc : expected.keySet() ) {
            Byte type = expected.get(loc);
            switch( type ) {
                case 'I':
                    if( !calc.get(loc).isSimpleInsertion() ) { return false; }
                    break;
                case 'D':
                    if( !calc.get(loc).isSimpleDeletion() ) { return false; }
                    break;
                case 'M':
                    if( !(calc.get(loc).isMNP() || calc.get(loc).isSNP()) ) { return false; }
                    break;
                default:
                    return false;
            }
        }
        return true;
    }

    @Test
    public void testReduceNumberOfAlternativeAllelesBasedOnHaplotypesScores(){

        // first have a list of alleles, one ref, several alt
        final Allele ref = Allele.create("A", true);
        final Allele altC = Allele.create("C", false);
        final Allele altT = Allele.create("T", false);
        final Allele altT2 = Allele.create("TT", false);
        final Allele altG = Allele.create("G", false);
        final Allele altT3 = Allele.create("TTT", false);

        // then create several haplotypes, assign ad-hoc scores
        final Haplotype hapRef = new Haplotype("AAAAA".getBytes());
        hapRef.setScore(Double.MAX_VALUE);

        // test case when both same best score and second best score are the same
        final Haplotype hapT = new Haplotype("TAAAA".getBytes());
        hapT.setScore(-2.0);
        final Haplotype hapTAnother = new Haplotype("TAAAT".getBytes());
        hapTAnother.setScore(-3.0);
        final Haplotype hapT2 = new Haplotype("TTAAA".getBytes());
        hapT2.setScore(-2.0);
        final Haplotype hapT2Another = new Haplotype("TTAAT".getBytes());
        hapT2Another.setScore(-3.0);

        final Haplotype hapC = new Haplotype("CAAAA".getBytes());
        hapC.setScore(-3.0);

        // for case when there's tie in highest haplotype score
        final Haplotype hapG = new Haplotype("GAAAA".getBytes());
        hapG.setScore(-3.0);
        final Haplotype hapGAnother = new Haplotype("GAAAG".getBytes());
        hapGAnother.setScore(-5.0);

        final Map<Allele, List<Haplotype>> alleleMapper = new LinkedHashMap<>();
        alleleMapper.put(ref, Arrays.asList(hapRef));
        alleleMapper.put(altC, Arrays.asList(hapC));
        alleleMapper.put(altT, Arrays.asList(hapT, hapTAnother));
        alleleMapper.put(altT2, Arrays.asList(hapT2, hapT2Another));
        alleleMapper.put(altG, Arrays.asList(hapG, hapGAnother));

        List<Allele> allelesToKeep = HaplotypeCallerGenotypingEngine.whichAllelesToKeepBasedonHapScores(alleleMapper, 5);
        Assert.assertEquals(allelesToKeep.size(), 5);

        Iterator<Allele> it = allelesToKeep.iterator();
        Assert.assertEquals(it.next(), ref);
        Assert.assertEquals(it.next(), altC);
        Assert.assertEquals(it.next(), altT);
        Assert.assertEquals(it.next(), altT2);
        Assert.assertEquals(it.next(), altG);

        allelesToKeep = HaplotypeCallerGenotypingEngine.whichAllelesToKeepBasedonHapScores(alleleMapper, 4);
        Assert.assertEquals(allelesToKeep.size(), 4);
        it = allelesToKeep.iterator();
        Assert.assertEquals(it.next(), ref);
        Assert.assertEquals(it.next(), altT);
        Assert.assertEquals(it.next(), altT2);
        Assert.assertEquals(it.next(), altG);

        allelesToKeep = HaplotypeCallerGenotypingEngine.whichAllelesToKeepBasedonHapScores(alleleMapper, 3);
        Assert.assertEquals(allelesToKeep.size(), 3);
        it = allelesToKeep.iterator();
        Assert.assertEquals(it.next(), ref);
        Assert.assertEquals(it.next(), altT);
        Assert.assertEquals(it.next(), altT2);

        allelesToKeep = HaplotypeCallerGenotypingEngine.whichAllelesToKeepBasedonHapScores(alleleMapper, 2);
        Assert.assertEquals(allelesToKeep.size(), 2);
        it = allelesToKeep.iterator();
        Assert.assertEquals(it.next(), ref);
        Assert.assertEquals(it.next(), altT);

        allelesToKeep = HaplotypeCallerGenotypingEngine.whichAllelesToKeepBasedonHapScores(alleleMapper, 1);
        Assert.assertEquals(allelesToKeep.size(), 1);
        it = allelesToKeep.iterator();
        Assert.assertEquals(it.next(), ref);

        // in the case of GGA mode there could be an allele with no haplotype support
        alleleMapper.put(altT3, new ArrayList<>());
        allelesToKeep = HaplotypeCallerGenotypingEngine.whichAllelesToKeepBasedonHapScores(alleleMapper, 1);
        Assert.assertEquals(allelesToKeep.size(), 1);
        it = allelesToKeep.iterator();
        Assert.assertEquals(it.next(), ref);
    }

    @Test
    public void testRemoveExcessiveAltAlleleFromVC(){
        final VariantContext originalVC = new VariantContextBuilder("source", "1", 1000000, 1000000, Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("C", false), Allele.create("G", false))).make();

        final VariantContext reducedVC = HaplotypeCallerGenotypingEngine.removeExcessAltAllelesFromVC(originalVC, Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("C", false)));

        Assert.assertEquals(reducedVC.getNAlleles(), 3);
        Assert.assertTrue(reducedVC.getAlleles().containsAll(Arrays.asList(Allele.create("A", true), Allele.create("T", false), Allele.create("C", false))));
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testMakeAnnotatedCallTrimmingAlleles(){
        List<Allele> alleles = Arrays.asList(Allele.create("AGGGGGGGGG", true), Allele.create("TGGGGGGGGG", false));
        List<Allele> mergedAlleles = Arrays.asList(Allele.create("AGGGGGGGGG", true), Allele.create("TGGGGGGGGG", false), Allele.create("A", false));
        AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<GATKRead, Allele>(SampleList.EMPTY_LIST, new IndexedAlleleList<Allele>(alleles), new HashMap<>());

        // Both a deletion and SNPs are present at this site
        final VariantContext originalVC = new VariantContextBuilder("source", "1", 1000000, 1000009, alleles).make();

        final List<FeatureInput<VariantContext>> features = Collections.emptyList();

        VariantContext reducedVC = HaplotypeCallerGenotypingEngine.makeAnnotatedCall("AGGGGGGGGG".getBytes(),
                new SimpleInterval(originalVC),
                new FeatureContext(),
                ArtificialReadUtils.createArtificialSamHeader(),
                originalVC,
                mergedAlleles.size(),
                likelihoods,
                originalVC,
                new VariantAnnotatorEngine(Collections.emptyList(), null, features, false, true), null);

        // Asserting that the two alleles were trimmed after calling removeExcessAltAlleles
        Assert.assertEquals(reducedVC.getNAlleles(), 2);
        Assert.assertTrue(reducedVC.getAlleles().containsAll(Arrays.asList(Allele.create("A", true), Allele.create("T", false))));
    }

    @Test
    public void testReplaceWithSpanDelVC() {
        final VariantContext snp = new VariantContextBuilder("source", "1", 1000000, 1000000,
                Arrays.asList(Allele.create("A", true),
                        Allele.create("G", false))).make();
        final VariantContext snpReplacement =
                HaplotypeCallerGenotypingEngine.replaceWithSpanDelVC(snp, Allele.create("A", true), 1000000);

        Assert.assertEquals(snpReplacement, snp);

        final VariantContext spanDel = new VariantContextBuilder("source", "1", 999995, 1000005,
                Arrays.asList(Allele.create("AAAAAAAAAAA", true),
                        Allele.create("A", false))).make();

        final VariantContext spanDelReplacement =
                HaplotypeCallerGenotypingEngine.replaceWithSpanDelVC(spanDel, Allele.create("A", true), 1000000);

        final VariantContext expectedSpanDelReplacement = new VariantContextBuilder("source", "1", 1000000, 1000000,
                Arrays.asList(Allele.create("A", true),
                        Allele.SPAN_DEL)).make();

        VariantContextTestUtils.assertVariantContextsAreEqual(spanDelReplacement,
                expectedSpanDelReplacement, Collections.emptyList(), Collections.emptyList());

        final VariantContext spanDelWithGt = new VariantContextBuilder("source", "1", 999995, 1000005,
                Arrays.asList(Allele.create("AAAAAAAAAAA", true),
                        Allele.create("A", false)))
                .genotypes(new GenotypeBuilder("S1",
                        Arrays.asList(spanDel.getAlleles().get(0), spanDel.getAlleles().get(1))).make()).make();

        final VariantContext spanDelWithGTReplacement =
                HaplotypeCallerGenotypingEngine.replaceWithSpanDelVC(spanDelWithGt, Allele.create("A", true), 1000000);

        VariantContextTestUtils.assertVariantContextsAreEqual(spanDelWithGTReplacement,
                expectedSpanDelReplacement, Collections.emptyList(), Collections.emptyList());

    }

    @DataProvider(name = "HomopolymerIndelTestProvider")
    public Object[][] makeHomopolymerIndelTests() {
        return new Object[][] {
                // ref sequence, loc in ref (0-based), ref allele, alt allele, hpolThreshold, expected result
                // Homopolymer insertion in long A run (10 As) - should be eligible with threshold 8
                { "GGGAAAAAAAAAAGGG", 3, "A", "AA", 8, true },
                // Homopolymer deletion in long A run (10 As) - should be eligible with threshold 8
                { "GGGAAAAAAAAAAGGG", 3, "AA", "A", 8, true },
                // Homopolymer insertion in short A run (5 As) - should NOT be eligible with threshold 8
                { "GGGAAAAAGGG", 3, "A", "AA", 8, false },
                // Homopolymer deletion in short A run (5 As) - should NOT be eligible with threshold 8
                { "GGGAAAAAGGG", 3, "AA", "A", 8, false },
                // Homopolymer insertion in T run (12 Ts) - should be eligible with threshold 10
                { "CCCTTTTTTTTTTTTCCC", 3, "T", "TT", 10, true },
                // Homopolymer deletion in T run (12 Ts) - should be eligible with threshold 10
                { "CCCTTTTTTTTTTTTCCC", 3, "TT", "T", 10, true },
                // SNP in homopolymer region - should NOT be eligible (not an indel)
                { "GGGAAAAAAAAAAGGG", 3, "A", "T", 8, false },
                // Indel not matching homopolymer base - should NOT be eligible
                { "GGGAAAAAAAAAAGGG", 3, "A", "AT", 8, false },
                // Homopolymer insertion with exactly threshold length - should be eligible
                { "GGGAAAAAAAAGGG", 3, "A", "AA", 8, true },
                // Homopolymer insertion with one less than threshold length - should NOT be eligible
                { "GGGAAAAAAAGGG", 3, "A", "AA", 8, false },
                // Null VC - should return false
                { "GGGAAAAAAAAAAGGG", 3, null, null, 8, false },
        };
    }

    @Test(dataProvider = "HomopolymerIndelTestProvider")
    public void testIsEligibleHomopolymerIndel(final String refSequence, final int loc,
                                               final String refAllele, final String altAllele,
                                               final int hpolThreshold, final boolean expected) {
        final byte[] ref = refSequence.getBytes();
        // Create DragstrReferenceAnalyzer with period 1 for homopolymer analysis
        final DragstrReferenceAnalyzer dragstrs = DragstrReferenceAnalyzer.of(ref, 0, ref.length, 1);

        final VariantContext vc;
        if (refAllele == null || altAllele == null) {
            vc = null;
        } else {
            vc = new VariantContextBuilder("source", "1", 1000000, 1000000 + refAllele.length() - 1,
                    Arrays.asList(Allele.create(refAllele, true), Allele.create(altAllele, false))).make();
        }

        final boolean result = GenotypingEngine.isEligibleHomopolymerIndel(vc, loc, dragstrs, hpolThreshold);
        Assert.assertEquals(result, expected,
                String.format("Failed for ref=%s, loc=%d, refAllele=%s, altAllele=%s, threshold=%d",
                        refSequence, loc, refAllele, altAllele, hpolThreshold));
    }

    @Test
    public void testIsHmerIndelViaIsEligibleHomopolymerIndel() {
        // Test homopolymer indel allele detection indirectly through isEligibleHomopolymerIndel
        // A 10-A homopolymer run
        final byte[] hpolRef = "GGGAAAAAAAAAAGGG".getBytes();
        final DragstrReferenceAnalyzer dragstrs = DragstrReferenceAnalyzer.of(hpolRef, 0, hpolRef.length, 1);

        // Position 3 is the first A in the homopolymer run
        final int loc = 3;

        // Test homopolymer A insertion (AAA allele) - should match A base and be eligible
        final VariantContext hmerInsertionA = new VariantContextBuilder("source", "1", 1000000, 1000000,
                Arrays.asList(Allele.create("A", true), Allele.create("AAA", false))).make();
        Assert.assertTrue(GenotypingEngine.isEligibleHomopolymerIndel(hmerInsertionA, loc, dragstrs, 8),
                "Homopolymer A insertion should be eligible in A homopolymer run");

        // Test mixed allele insertion (ATG allele) - should NOT be eligible (not a homopolymer allele)
        final VariantContext mixedInsertion = new VariantContextBuilder("source", "1", 1000000, 1000000,
                Arrays.asList(Allele.create("A", true), Allele.create("ATG", false))).make();
        Assert.assertFalse(GenotypingEngine.isEligibleHomopolymerIndel(mixedInsertion, loc, dragstrs, 8),
                "Mixed allele insertion should NOT be eligible");

        // Test homopolymer T insertion in A run - should NOT be eligible (wrong base)
        final VariantContext hmerInsertionT = new VariantContextBuilder("source", "1", 1000000, 1000000,
                Arrays.asList(Allele.create("A", true), Allele.create("TT", false))).make();
        Assert.assertFalse(GenotypingEngine.isEligibleHomopolymerIndel(hmerInsertionT, loc, dragstrs, 8),
                "Homopolymer T insertion should NOT be eligible in A homopolymer run");
    }

    @Test
    public void testHomopolymerGenotypingThresholdNonZero() {
        // Test that homopolymer indels are correctly identified with various thresholds

        // 15 As in a row - should be eligible with reasonable thresholds
        final byte[] longHpolRef = "CCCAAAAAAAAAAAAAAACCC".getBytes();
        final DragstrReferenceAnalyzer dragstrsLong = DragstrReferenceAnalyzer.of(longHpolRef, 0, longHpolRef.length, 1);

        // Verify the period and repeat length at position 3 (first A)
        Assert.assertEquals(dragstrsLong.period(3), 1, "Period should be 1 for homopolymer");
        Assert.assertTrue(dragstrsLong.repeatLength(3) >= 10, "Repeat length should be >= 10 for long homopolymer");

        // Create an insertion in the homopolymer
        final VariantContext hpolInsertion = new VariantContextBuilder("source", "1", 1000000, 1000000,
                Arrays.asList(Allele.create("A", true), Allele.create("AA", false))).make();

        // Should be eligible with threshold of 8
        Assert.assertTrue(GenotypingEngine.isEligibleHomopolymerIndel(hpolInsertion, 3, dragstrsLong, 8));
        // Should be eligible with threshold of 10
        Assert.assertTrue(GenotypingEngine.isEligibleHomopolymerIndel(hpolInsertion, 3, dragstrsLong, 10));
        // Should NOT be eligible with threshold of 20 (exceeds homopolymer length)
        Assert.assertFalse(GenotypingEngine.isEligibleHomopolymerIndel(hpolInsertion, 3, dragstrsLong, 20));

        // Create a deletion in the homopolymer
        final VariantContext hpolDeletion = new VariantContextBuilder("source", "1", 1000000, 1000001,
                Arrays.asList(Allele.create("AA", true), Allele.create("A", false))).make();

        // Should be eligible with threshold of 8
        Assert.assertTrue(GenotypingEngine.isEligibleHomopolymerIndel(hpolDeletion, 3, dragstrsLong, 8));

        // 5 As in a row - should NOT be eligible with threshold of 8
        final byte[] shortHpolRef = "CCCAAAAACCC".getBytes();
        final DragstrReferenceAnalyzer dragstrsShort = DragstrReferenceAnalyzer.of(shortHpolRef, 0, shortHpolRef.length, 1);

        Assert.assertFalse(GenotypingEngine.isEligibleHomopolymerIndel(hpolInsertion, 3, dragstrsShort, 8));
    }
}
