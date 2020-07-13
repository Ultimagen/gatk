package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.ultimagenomics.flow_based_read.read.FlowBasedHaplotype;

import java.util.*;
import java.util.stream.Collectors;

//class implementing a boolean matrix with rows - haplotypes, columns, LocationAlleles and T/F with coccurence
public class HaplotypeAlleleMatrix {
    List<Haplotype> haplotypeList;

    boolean[][] coocurrenceMatrix;
    Map<LocationAndAlleles, Integer> variant2Col;
    LocationAndAlleles[] col2Variant;

    final int n_haplotypes;
    final int n_variants;


    public HaplotypeAlleleMatrix(List<Haplotype> _haplotypeList) {
        haplotypeList = _haplotypeList;
        variant2Col = new HashMap<>();
        int col_count = 0;

        for (Haplotype h : haplotypeList) {
            Collection<VariantContext> variants = h.getEventMap().values();
            for (VariantContext var : variants) {
                if (!variant2Col.containsKey(new LocationAndAlleles(var.getStart(), var.getAlleles()))) {
                    variant2Col.put(new LocationAndAlleles(var.getStart(), var.getAlleles()), col_count);
                    col_count++;
                }
            }
        }

        col2Variant = new LocationAndAlleles[variant2Col.size()];

        for (LocationAndAlleles k : variant2Col.keySet()) {
            col2Variant[variant2Col.get(k)] =  k;
        }

        coocurrenceMatrix = new boolean[haplotypeList.size()][col_count];

        for (int i = 0; i < haplotypeList.size(); i++) {
            Collection<VariantContext> variants = haplotypeList.get(i).getEventMap().values();
            for (VariantContext var : variants) {
                coocurrenceMatrix[i][variant2Col.get(new LocationAndAlleles(var.getStart(), var.getAlleles()))] = true;
            }
        }
        n_haplotypes = haplotypeList.size();
        n_variants = col_count;
    }


    public List<Pair<LocationAndAlleles, LocationAndAlleles>> nonCoOcurringVariants() {

        List<Pair<Integer, Integer>> result = new ArrayList<>();
        for (int i = 0; i < n_variants; i++) {
            for (int j = i + 1; j < n_variants; j++) {
                boolean flag = false;
                for (int r = 0; r < n_haplotypes; r++) {
                    if (coocurrenceMatrix[r][i] & coocurrenceMatrix[r][j]) {
                        flag = true;
                        break;
                    }
                }
                if (!flag) {
                    result.add(ImmutablePair.of(i,j));
                }

            }
        }
        List<Pair<LocationAndAlleles,LocationAndAlleles>> vc_result = new ArrayList<>();
        for (Pair<Integer, Integer> res: result) {
            vc_result.add(ImmutablePair.of(col2Variant[res.getLeft()], col2Variant[res.getRight()]));
        }
        return vc_result;
    }


    static public List<Pair<LocationAndAlleles, LocationAndAlleles>>
            filterExclusivePairsByDistance(List<Pair<LocationAndAlleles, LocationAndAlleles>> allelePairs,
                                           int max_dist) {
        allelePairs.removeIf(v -> Math.abs(v.getLeft().getLoc() - v.getRight().getLoc())>max_dist);
        allelePairs.removeIf(v -> Math.abs(v.getLeft().getLoc() - v.getRight().getLoc())==0);

        return allelePairs;
    }

    public List<Pair<LocationAndAlleles, LocationAndAlleles>>
            filterSameUpToHmerPairs(List<Pair<LocationAndAlleles, LocationAndAlleles>> allelePairs, int activeWindowStart) {
        List<Pair<LocationAndAlleles, LocationAndAlleles>> result = new ArrayList<>();
        List<Haplotype> modifiedHaplotypes = insertAllelesIntoReference(activeWindowStart);
        List<FlowBasedHaplotype> fbh = haplotypes2FlowHaplotypes(modifiedHaplotypes);
        Set<Pair<Integer, Integer>> upToHmerList = upToHmer(fbh);
        for (Pair<LocationAndAlleles, LocationAndAlleles> p : allelePairs) {
            int first = variant2Col.get(p.getLeft());
            int second = variant2Col.get(p.getRight());
            if (upToHmerList.contains(ImmutablePair.of(first,second)))
            {
                result.add(p);
            }
        }
        return result;
    }

    static public Map<LocationAndAlleles, Set<LocationAndAlleles>> getExclusiveAlleleMap(List<Pair<LocationAndAlleles, LocationAndAlleles>> input) {
        Map<LocationAndAlleles, Set<LocationAndAlleles>> result = new LinkedHashMap<>();
        for ( Pair<LocationAndAlleles, LocationAndAlleles> v: input ) {
            LocationAndAlleles key = v.getLeft();
            LocationAndAlleles value = v.getRight();
            if (result.containsKey(key)) {
                result.get(key).add(value);
            } else {
                result.put(key, new HashSet<>());
                result.get(key).add(value);
            }

            key = v.getRight();
            value = v.getLeft();
            if (result.containsKey(key)) {
                result.get(key).add(value);
            } else {
                result.put(key, new HashSet<>());
                result.get(key).add(value);
            }
        }
        return result;
    }

    private Haplotype findReferenceHaplotype() {
        for (Haplotype h: haplotypeList ) {
            if (h.isReference()) {
                return h;
            }
        }
        return null;
    }


    private List<Haplotype> insertAllelesIntoReference(int activeWindowStart) {
        final Haplotype reference = findReferenceHaplotype();
        List<Haplotype> modifiedHaplotypes = new ArrayList<>();

        for (int i = 0 ; i < col2Variant.length; i++) {
            LocationAndAlleles tmp = col2Variant[i];
            int commonPrefixShift = commonPrefixLength(tmp);

            List<Allele> allelesList = tmp.getAlleles();
            if (allelesList.size()>2) {
                throw new RuntimeException("More than two alleles extracted from haplotype");
            }
            Allele ref=null;
            Allele alt=null;

            for ( Allele al: allelesList ) {
                if (al.isNonReference()){
                    alt = al;
                } else if (al.isReference()) {
                    ref = al;
                }
            }

            Haplotype tmpHap = reference.insertAllele(ref, alt,
                    tmp.getLoc() - activeWindowStart,
                    tmp.getLoc() - activeWindowStart,
                    commonPrefixShift);
            modifiedHaplotypes.add(tmpHap);

        }
        return modifiedHaplotypes;
    }

    //removes common prefix from indels
    private int commonPrefixLength(LocationAndAlleles input){
        List<Allele> alleleList = input.getAlleles();

        for (Allele al: alleleList){
            if (al.isSymbolic())
                return 0;
        }
        if (alleleList.size() <=1)
            return 0;

        int minLen = alleleList.get(0).length();
        for (int i = 1; i < alleleList.size(); i++){
            if (alleleList.get(i).length() < minLen)
                minLen = alleleList.get(i).length();
        }

        boolean foundNonCommonBase = false;
        int shift = 0 ;
        for (int i = input.getLoc(); i < input.getLoc()+minLen; i++){
            byte base = alleleList.get(0).getBases()[i-input.getLoc()];
            for (int j = 1; j < alleleList.size(); j++){
                if (alleleList.get(j).getBases()[i-input.getLoc()] != base){
                    foundNonCommonBase = true;
                }
            }
            if (foundNonCommonBase) {
                shift = i - input.getLoc();
                break;
            }
        }
        if (!foundNonCommonBase) {
            shift = minLen;
        }

        return shift;

    }
    private List<FlowBasedHaplotype> haplotypes2FlowHaplotypes(List<Haplotype> haps) {
        List<FlowBasedHaplotype> flowBasedHaplotypes = haps.stream().map(h -> new FlowBasedHaplotype(h, "GCAT")).collect(Collectors.toList());

        return flowBasedHaplotypes;
    }

    private Set<Pair<Integer, Integer>> upToHmer(List<FlowBasedHaplotype> fbhaps) {
        Set<Pair<Integer,Integer>> result = new HashSet<>();
        for (int i = 0 ; i < fbhaps.size(); i++ ) {
            for (int j = i + 1; j < fbhaps.size(); j++) {
                if (fbhaps.get(i).equalUpToHmerChange(fbhaps.get(j))) {
                    result.add(ImmutablePair.of(i, j));
                    result.add(ImmutablePair.of(j, i));
                }
            }
        }
        return result;
    }


}

