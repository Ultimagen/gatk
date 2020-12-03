package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.AnnotationUtils;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;

@DocumentedFeature(groupName= HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Allele-specific length of homopolymer (longer of either ref or alt)")
public class AS_HmerLength extends InfoFieldAnnotation implements AS_StandardAnnotation, AlleleSpecificAnnotation, ReducibleAnnotation {
    @Override
    public Map<String, Object> annotate(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        return Collections.emptyMap();
        //TODO: make sure regular HC still outputs the annotation
    }

    @Override
    public String getPrimaryRawKey() {
        return GATKVCFConstants.AS_RAW_HMER_LENGTH_KEY;
    }

    @Override
    public boolean hasSecondaryRawKeys() {
        return false;
    }

    @Override
    public List<String> getSecondaryRawKeys() {
        return null;
    }

    @Override
    public Map<String, Object> annotateRawData(ReferenceContext ref, VariantContext vc, AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(ref, "ref is null");
        if (!ref.hasBackingDataSource()) {
            throw new UserException("Reference is required to calculate AS_HmerLength annotation.");
        }

        final List<VariantContext> splitVcs = GATKVariantContextUtils.splitVariantContextToBiallelics(vc, true, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL, false);
        final List<String> hmerLengths = new LinkedList<>();

        for (final VariantContext splitVariant : splitVcs) {
            final Allele altAllele = splitVariant.getAlternateAllele(0);

            if (splitVariant.isSNP()) {
                hmerLengths.add(String.valueOf(1));
            } else if (splitVariant.isMNP()) {
                hmerLengths.add(String.valueOf(0));
            } else if (splitVariant.isSymbolic()) {
                hmerLengths.add("");
            } else {
                final byte base = ref.getBases()[1];
                int altHmerLength = 0;
                for (int i=1; i < altAllele.getBases().length; i++) {
                    final byte curBase = altAllele.getBases()[i];
                    if (curBase == base) {
                        altHmerLength++;
                    } else {
                        altHmerLength = 0;
                    }
                }

                int refHmerLength = 0;
                boolean isFirstBaseOfRef = true;
                for (Byte refBase : ref) {
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

                int finalHmerLength = variantIsHmer(splitVariant) ? refHmerLength + altHmerLength : 0;
                hmerLengths.add(String.valueOf(finalHmerLength));
            }
        }

        final Map<String, Object> map = new HashMap<>();
        map.put(getPrimaryRawKey(), String.join(AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM, hmerLengths));
        return map;
    }

    @Override
    public Map<String, Object> combineRawData(List<Allele> allelesList, List<ReducibleAnnotationData<?>> listOfRawData) {
        Map<Allele, String> hmerLengths = new HashMap<>();

        for (final ReducibleAnnotationData currentValue : listOfRawData) {
            List<String> currentHmerLengths = AnnotationUtils.getAlleleLengthListOfString(currentValue.rawData);

            for (int i=0; i<allelesList.size(); i++) {
                Allele a = allelesList.get(i);
                if (!a.isReference()) {
                    if (hmerLengths.containsKey(a)) {
                        if (!currentValue.attributeMap.get(a).equals(hmerLengths.get(a))) {
                            throw new UserException("Hmer lengths for the same allele must match.");
                        }
                    } else {
                        //we don't have an entry for the ref allele, so subtract 1 from the index
                        hmerLengths.put(a, currentHmerLengths.get(i-1));
                    }
                }
            }
        }

        List<String> hmerLengthsInAlleleOrder = new LinkedList<>();
        for(Allele a : allelesList) {
            if (!a.isReference()) {
                hmerLengthsInAlleleOrder.add(hmerLengths.get(a));
            }
        }
        return Collections.singletonMap(getPrimaryRawKey(), String.join(AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM, hmerLengthsInAlleleOrder));
        //return null;
    }

    /*
    // To "finalize" this annotation, just replace the raw delimiter in the list to the reduced delimiter
     */
    @Override
    public Map<String, Object> finalizeRawData(VariantContext vc, VariantContext originalVC) {
        if (!vc.hasAttribute(getPrimaryRawKey())) {
            return new HashMap<>();
        }
        String rawHmerLengthData = vc.getAttributeAsString(getPrimaryRawKey(),null);
        if (rawHmerLengthData == null) {
            return new HashMap<>();
        }
        AlleleSpecificAnnotationData<List<Integer>> myData = new AlleleSpecificAnnotationData<>(originalVC.getAlleles(), rawHmerLengthData);

        Map<String, Object> returnMap = new HashMap<>();
        returnMap.put(getKeyNames().get(0), myData.rawData.replace(AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM, AnnotationUtils.ALLELE_SPECIFIC_REDUCED_DELIM));
        returnMap.put(getPrimaryRawKey(), myData.rawData);  //this is in case raw annotations are requested
        return returnMap;
    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.AS_HMER_LENGTH_KEY);
    }

    // Checks if biallelic INDEL is a homopolymer (need to check both ref and alt)
    private Boolean variantIsHmer(VariantContext vc) {
        Allele refAllele = vc.getReference();
        Allele altAllele = vc.getAlternateAllele(0);
        Allele longerAllele = refAllele.length() > altAllele.length() ? refAllele : altAllele;
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

}
