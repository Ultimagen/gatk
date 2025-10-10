package org.broadinstitute.hellbender.tools.walkers.annotator.flow;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.FlowBasedAlignmentArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Unit tests for HmerIndelBiasAnnotation
 *
 * This test class validates the basic functionality of the HmerIndelBiasAnnotation class,
 * focusing on interface compliance and basic properties that can be tested without
 * complex mocking of flow-based reads and reference contexts.
 */
public class HmerIndelBiasTest extends GATKBaseTest {
    private static final String SAMPLE = "NA12878";
    public static final String TEST_FILES_DIR = toolsTestDir + "walkers/annotator/flow/";

    @Test
    public void testDescription(){
        String[] constants = {GATKVCFConstants.FLOW_HMER_INDEL_BIAS_KEY};
        VCFFormatHeaderLine[] hlines    = {GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.FLOW_HMER_INDEL_BIAS_KEY)};
        Assert.assertEquals(new HmerIndelBias().getKeyNames(), new ArrayList<>(Arrays.asList(constants)));
        Assert.assertEquals(new HmerIndelBias().getDescriptions(), new ArrayList<>(Arrays.asList(hlines)));
    }

    @DataProvider(name = "hmerIndelBiasCases")
    public Object[][] hmerIndelBiasCases() {
        return new Object[][] {
                {
                        Allele.create("C", true), // REF
                        Allele.create("G"),       // ALT
                        TEST_FILES_DIR + "no_variant.bam", // BAM
                        new SimpleInterval("ref", 150, 150), // Location
                        false,
                        null,// Expected output
                        "NoVariant"
                },
                {
                        Allele.create("G", true),
                        Allele.create("GT"),
                        TEST_FILES_DIR + "hmer_ins_122.bam",
                        new SimpleInterval("ref", 122, 122),
                        true,
                        "8,16|0,0",
                        "HmerInsVariant"
                },
                {
                        Allele.create("GT", true),
                        Allele.create("G"),
                        TEST_FILES_DIR + "hmer_del_122.bam",
                        new SimpleInterval("ref", 122, 123),
                        true,
                        "10,11|0,0",
                        "HmerDelVariant"
                },

                {
                        Allele.create("C", true),
                        Allele.create("G"),
                        TEST_FILES_DIR + "snp_150.bam",
                        new SimpleInterval("ref", 150, 150),
                        false,
                        null,
                        "SNPVariant"
                },
                {
                        Allele.create("G", true),
                        Allele.create("GA"),
                        TEST_FILES_DIR + "indel_150.bam",
                        new SimpleInterval("ref", 150, 150),
                        false,
                        null,
                        "NonHmerVariant"
                }


        };
    }

    @Test(dataProvider = "hmerIndelBiasCases")
    public void testHmerIndelBiasAnnotation(
            Allele REF,
            Allele ALT,
            String bamPath,
            SimpleInterval interval,
            boolean outputExpected,
            String expectedOutput,
            String testName
    ) {
        final List<Allele> ALLELES = Arrays.asList(REF, ALT);

        final CachingIndexedFastaSequenceFile referenceReader =
                new CachingIndexedFastaSequenceFile(IOUtils.getPath(TEST_FILES_DIR + "ref.fa"));

        final ReadsDataSource readsDataSource =
                new ReadsPathDataSource(IOUtils.getPath(bamPath));

        List<GATKRead> refReads = new ArrayList<>();
        List<GATKRead> altReads = new ArrayList<>();
        for (GATKRead read : readsDataSource) {
            refReads.add(new FlowBasedRead(read, "TGCA", 20, new FlowBasedAlignmentArgumentCollection()));
        }

        HmerIndelBias annotator = new HmerIndelBias();

        final int dpDepth = 50;
        final Genotype gAC = new GenotypeBuilder(SAMPLE, ALLELES).DP(dpDepth).make();
        final double log10PError = -5;

        final AlleleLikelihoods<GATKRead, Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(SAMPLE, refReads, altReads, -100.0, -100.0, REF, ALT);

        final VariantContext vc = new VariantContextBuilder("test", "ref", interval.getStart(), interval.getEnd(), ALLELES)
                .log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();

        final GenotypeBuilder gb = new GenotypeBuilder(gAC);
        SimpleInterval referenceInterval = new SimpleInterval(vc.getContig(), vc.getStart()-20, vc.getEnd()+20);
        ReferenceContext refContext = new ReferenceContext(new ReferenceFileSource(IOUtils.getPath(TEST_FILES_DIR + "ref.fa")), referenceInterval);

        annotator.annotate(refContext, vc, gAC, gb, likelihoods);
        if (outputExpected){
            Assert.assertEquals(gb.make().hasExtendedAttribute("X_HIB"), true, "Test: " + testName);
            Assert.assertEquals(gb.make().getExtendedAttribute("X_HIB"), expectedOutput, "Test: " + testName);

        } else {
            Assert.assertEquals(gb.make().hasExtendedAttribute("X_HIB"), false, "Test: " + testName);
        }

        referenceReader.close();
        readsDataSource.close();
    }

}
