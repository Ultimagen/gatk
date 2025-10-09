package org.broadinstitute.hellbender.tools.walkers.annotator.flow;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.testutils.ArtificialAnnotationUtils;
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
    private static final Allele REF = Allele.create("C", true);
    private static final Allele ALT = Allele.create("G");
    private static final Allele REF2 = Allele.create("T", true);
    private static final Allele ALT2 = Allele.create("TC");
    private static final List<Allele> ALLELES = Arrays.asList(REF, ALT);
    private static final List<Allele> ALLELES2 = Arrays.asList(REF2, ALT2);
    private static final String SAMPLE = "NA12878";
    public static final String TEST_FILES_DIR = toolsTestDir + "walkers/annotator/flow/";

    @Test
    public void testDescription(){
        String[] constants = {GATKVCFConstants.FLOW_HMER_INDEL_BIAS_KEY};
        VCFFormatHeaderLine[] hlines    = {GATKVCFHeaderLines.getFormatLine(GATKVCFConstants.FLOW_HMER_INDEL_BIAS_KEY)};
        Assert.assertEquals(new HmerIndelBias().getKeyNames(), new ArrayList<>(Arrays.asList(constants)));
        Assert.assertEquals(new HmerIndelBias().getDescriptions(), new ArrayList<>(Arrays.asList(hlines)));
    }

    @Test
    public void testAnnotationNoVariant(){
        // Read reference genome from exampleFASTA.fasta

        final CachingIndexedFastaSequenceFile referenceReader = 
            new CachingIndexedFastaSequenceFile(IOUtils.getPath(TEST_FILES_DIR + "ref.fa"));
        
        // Read BAM file no_variant.bam from the flow test resources directory
        final ReadsDataSource readsDataSource = 
            new ReadsPathDataSource(IOUtils.getPath(TEST_FILES_DIR + "no_variant.bam"));
        
        // Convert reads to FlowBasedReads and separate into ref and alt reads
        List<GATKRead> refReads = new ArrayList<>();
        List<GATKRead> altReads = new ArrayList<>();
        for (GATKRead read : readsDataSource) {
            if (read instanceof FlowBasedRead) {
                // For this test, we'll put all reads in refReads since we don't have variant calling logic
                refReads.add(read);
            }
        }
        
        // Test the annotation with the created reads
        HmerIndelBias annotator = new HmerIndelBias();

        final int dpDepth = 50; //Note: using a different value on purpose so that we can check that reads are preferred over DP
        final Genotype gAC = new GenotypeBuilder(SAMPLE, ALLELES).DP(dpDepth).make();
        final double log10PError = -5;

        final AlleleLikelihoods<GATKRead, Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(SAMPLE, refReads, altReads, -100.0, -100.0, REF, ALT);

        final VariantContext vc = new VariantContextBuilder("test", "ref", 150, 150, ALLELES).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();

        final GenotypeBuilder gb = new GenotypeBuilder(gAC);
        final SimpleInterval interval = new SimpleInterval("ref", 140, 160);
        ReferenceContext refContext = new ReferenceContext(new ReferenceFileSource(IOUtils.getPath(TEST_FILES_DIR + "ref.fa")), interval);

        annotator.annotate(refContext, vc, gAC, gb, likelihoods);

        Assert.assertEquals(gb.make().hasExtendedAttribute("X_HIB"),false);

        // Close resources
        referenceReader.close();
        readsDataSource.close();
    }

    @Test
    public void testAnnotationNoHmerVariant(){
        // Read reference genome from exampleFASTA.fasta

        final CachingIndexedFastaSequenceFile referenceReader =
                new CachingIndexedFastaSequenceFile(IOUtils.getPath(TEST_FILES_DIR + "ref.fa"));

        // Read BAM file no_variant.bam from the flow test resources directory
        final ReadsDataSource readsDataSource =
                new ReadsPathDataSource(IOUtils.getPath(TEST_FILES_DIR + "no_variant.bam"));

        // Convert reads to FlowBasedReads and separate into ref and alt reads
        List<GATKRead> refReads = new ArrayList<>();
        List<GATKRead> altReads = new ArrayList<>();
        for (GATKRead read : readsDataSource) {
            if (read instanceof FlowBasedRead) {
                // For this test, we'll put all reads in refReads since we don't have variant calling logic
                refReads.add(read);
            }
        }

        // Test the annotation with the created reads
        HmerIndelBias annotator = new HmerIndelBias();

        final int dpDepth = 50; //Note: using a different value on purpose so that we can check that reads are preferred over DP
        final Genotype gAC = new GenotypeBuilder(SAMPLE, ALLELES2).DP(dpDepth).make();
        final double log10PError = -5;

        final AlleleLikelihoods<GATKRead, Allele> likelihoods =
                ArtificialAnnotationUtils.makeLikelihoods(SAMPLE, refReads, altReads, -100.0, -100.0, REF, ALT);

        final VariantContext vc = new VariantContextBuilder("test", "ref", 10, 10, ALLELES2).log10PError(log10PError).genotypes(Arrays.asList(gAC)).make();

        final GenotypeBuilder gb = new GenotypeBuilder(gAC);
        final SimpleInterval interval = new SimpleInterval("ref", 1, 20);
        ReferenceContext refContext = new ReferenceContext(new ReferenceFileSource(IOUtils.getPath(TEST_FILES_DIR + "ref.fa")), interval);

        annotator.annotate(refContext, vc, gAC, gb, likelihoods);

        Assert.assertEquals(gb.make().hasExtendedAttribute("X_HIB"),false);

        // Close resources
        referenceReader.close();
        readsDataSource.close();
    }

}
