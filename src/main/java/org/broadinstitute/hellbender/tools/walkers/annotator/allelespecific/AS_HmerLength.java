package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.tools.walkers.variantutils.LeftAlignAndTrimVariants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

/**
 * Allele specific annotation of length of longest homopolymer. In the case of deletions, the longest homopolymer is
 * however long the reference h-mer is (before the deletion). For insertions it is the length of the h-mer in the reference
 * plus the length of the insertion. Indels with multiple h-mer changes are excluded (given a length of 0).
 *
 * This is not a reducible annotation, meaning it must be calculated after combining multiple GVCFS, or when generating a
 * single VCF. This is done by adding -A AS_HmerLength when calling HaplotypeCaller in VCF mode, or when calling
 * GenotypeGvcfs (GnarlyGenotyper will be adding -A arg at some point in the future, but doesn't work today). If the annotation
 * is generated in a GVCF, there is currently no way to combine it when JointCalling. If GenomicsDB is used, it will automatically
 * drop the annotation.
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Allele-specific length of homopolymer (longer of either ref or alt)")
public class AS_HmerLength extends InfoFieldAnnotation implements AS_StandardAnnotation, AlleleSpecificAnnotation {
    @Override
    public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(ref, "ref is null");
        if (!ref.hasBackingDataSource()) {
            throw new UserException("Reference is required to calculate AS_HmerLength annotation.");
        }

        final List<VariantContext> splitVcs = GATKVariantContextUtils.splitVariantContextToBiallelics(vc, true, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, false);
        final SimpleInterval refWindowFromVcStart = new SimpleInterval (vc.getContig(), vc.getStart(), ref.getWindow().getEnd());
        final List<String> hmerLengths = new LinkedList<>();

        for (final VariantContext splitVariant : splitVcs) {
            //splitVariant is now biallelic so we can get the first (only) alt allele.
            final Allele altAllele = splitVariant.getAlternateAllele(0);

            if (splitVariant.isSNP()) {
                hmerLengths.add(String.valueOf(0));
            } else if (splitVariant.isMNP()) {
                hmerLengths.add(String.valueOf(0));
            } else if (splitVariant.isSymbolic()) {
                hmerLengths.add("");
            } else {
                if (variantIsHmer(splitVariant) && !indelIsLeftAligned(splitVariant, ref)) {
                    throw new UserException(String.format("H-mer indel at %s is not left aligned.", vc.toString()));
                }
                //INDELs are anchored on the left by an unchanged base, so the h-mer base in question is the second base from variant start in the longer allele.
                final Allele longerAllele = splitVariant.getReference().length() > altAllele.length() ? splitVariant.getReference() : altAllele;
                final byte base = longerAllele.getBases()[1];
                int altHmerLength = 0;
                for (int i=1; i < altAllele.getBases().length; i++) {
                    final byte curBase = altAllele.getBases()[i];
                    if (curBase == base) {
                        altHmerLength++;
                    } else {
                        altHmerLength = 0;
                        break;
                    }
                }

                int refHmerLength = 0;
                boolean isFirstBaseOfRef = true;
                for (Byte refBase : ref.getBases(refWindowFromVcStart)) {
                    if (isFirstBaseOfRef) {
                        isFirstBaseOfRef = false;
                        continue;
                    }
                    if (base == refBase) {
                        refHmerLength++;
                    } else {
                        break;
                    }
                }

                //Still need to check if variant is actually an Hmer since we only checked the alt allele, not the
                // ref allele (just the reference itself).
                int finalHmerLength = variantIsHmer(splitVariant) ? refHmerLength + altHmerLength : 0;
                hmerLengths.add(String.valueOf(finalHmerLength));
            }
        }

        final Map<String, Object> map = new HashMap<>();
        map.put(getKeyNames().get(0), String.join(AnnotationUtils.ALLELE_SPECIFIC_REDUCED_DELIM, hmerLengths));
        return map;
    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.AS_HMER_LENGTH_KEY);
    }


    /**
     * Checks if biallelic INDEL is a homopolymer. Need to check either ref or alt are all one nucleotide (whichever is longer).
     *
     * @param vc Biallelic INDEL variant context. Must only have one alt allele.
     */
    private Boolean variantIsHmer(VariantContext vc) {
        Allele refAllele = vc.getReference();
        Allele altAllele = vc.getAlternateAllele(0);
        Allele longerAllele = refAllele.length() > altAllele.length() ? refAllele : altAllele;
        //INDELs are left aligned so start at the second base to test for an h-mer.
        byte base = longerAllele.getBases()[1];
        boolean isHmer = false;
        for (int i=1; i < longerAllele.getBases().length; i++) {
            final byte curBase = longerAllele.getBases()[i];
            if (curBase != base) {
                isHmer = false;
                break;
            } else {
                isHmer = true;
            }
        }
        return isHmer;
    }

    /**
     * Checks if bialleleic INDEL is left aligned
     *
     * @param vc Bialleliic INDEL variant context
     * @param ref reference context for this site
     * @return
     */
    private Boolean indelIsLeftAligned(VariantContext vc, ReferenceContext ref) {
        if (!vc.isBiallelic()){
            throw new GATKException(String.format("Expected VariantContext to be biallelic: %s", vc.toString()));
        }
        VariantContext leftAligned = GATKVariantContextUtils.leftAlignAndTrim(vc, ref, LeftAlignAndTrimVariants.DEFAULT_MAX_LEADING_BASES, false);
        return leftAligned.getReference() == vc.getReference() && leftAligned.getAlternateAllele(0) == leftAligned.getAlternateAllele(0);
    }

}
