package org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandBiasBySample;
import org.broadinstitute.hellbender.tools.walkers.annotator.StrandBiasTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StrandBiasTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AlleleSpecificAnnotationData;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotationData;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;


public abstract class AS_StrandBiasTestProbabilistic extends StrandBiasTest implements ReducibleAnnotation {
    private final static Logger logger = LogManager.getLogger(AS_StrandBiasTestProbabilistic.class);

    public static final int FORWARD = 0;
    public static final int REVERSE = 1;
    public static final String SPLIT_DELIM = "\\|"; //String.split takes a regex, so we need to escape the pipe
    public static final String PRINT_DELIM = "|";
    public static final String REDUCED_DELIM = ",";
    public static final int MIN_COUNT = 2;
    public static final double MIN_PVALUE = 1.0E-320;
    private final List<Double> ZERO_LIST = new ArrayList<>();

    @Override
    public String getPrimaryRawKey() { return GATKVCFConstants.AS_SBP_TABLE_KEY; }

    /**
     * @return true if annotation has secondary raw keys
     */
    @Override
    public boolean hasSecondaryRawKeys() {
        return false;
    }

    /**
     * Get additional raw key strings that are not the primary key
     *
     * @return may be null
     */
    @Override
    public List<String> getSecondaryRawKeys() {
        return null;
    }

    public AS_StrandBiasTestProbabilistic(){
        ZERO_LIST.add(0,0.0);
        ZERO_LIST.add(1,0.0);
    }

    protected String makeRawAnnotationString(final List<Allele> vcAlleles, final Map<Allele, List<Double>> perAlleleValues) {
        String annotationString = "";
        for (final Allele a : vcAlleles) {
            if (!annotationString.isEmpty()) {
                annotationString += PRINT_DELIM;
            }
            List<Double> alleleValues = perAlleleValues.get(a);
            if (alleleValues == null) {
                alleleValues = ZERO_LIST;
            }
            annotationString += encode(alleleValues);
        }
        return annotationString;
    }

    protected String encode(List<Double> alleleValues) {
        String annotationString = "";
        for (int j =0; j < alleleValues.size(); j++) {
            annotationString += alleleValues.get(j);
            if (j < alleleValues.size()-1) {
                annotationString += ",";
            }
        }
        return annotationString;
    }

    /**
     * Uses the ReadLikelihoods map to generate a 2x2 strand contingency table by counting the total read support for each
     * allele in either the forward or reverse direction.
     *
     * @param ref the reference context for this annotation
     * @param vc the variant context to annotate
     * @param likelihoods likelihoods indexed by sample, allele, and read within sample
     * @return
     */
    @Override
    public Map<String, Object> annotateRawData(final ReferenceContext ref,
                                               final VariantContext vc,
                                               final AlleleLikelihoods<GATKRead, Allele> likelihoods ) {

        //for allele-specific annotations we only call from HC and we only use likelihoods
        if ( likelihoods == null || !likelihoods.hasFilledLikelihoods()) {
            return Collections.emptyMap();
        }
        // calculate the annotation from the likelihoods
        // likelihoods can come from HaplotypeCaller call to VariantAnnotatorEngine
        final Map<String, Object> annotations = new HashMap<>();
        final ReducibleAnnotationData<List<Double>> myData = new AlleleSpecificAnnotationData<>(vc.getAlleles(),null);
        getStrandCountsFromLikelihoodMap(vc, likelihoods, myData, MIN_COUNT);
        final Map<Allele, List<Double>> perAlleleValues = myData.getAttributeMap();
        final String annotationString = makeRawAnnotationString(vc.getAlleles(), perAlleleValues);
        annotations.put(getPrimaryRawKey(), annotationString);
        return annotations;
    }


    public void getStrandCountsFromLikelihoodMap( final VariantContext vc,
                                                  final AlleleLikelihoods<GATKRead, Allele> likelihoods,
                                                  final ReducibleAnnotationData<List<Double>> perAlleleValues,
                                                  final int minCount) {
        if( likelihoods == null || vc == null ) {
            return;
        }

        final Allele ref = vc.getReference();
        final List<Allele> allAlts = vc.getAlternateAlleles();

        for (final String sample : likelihoods.samples()) {
            final ReducibleAnnotationData<List<Double>> sampleTable = new AlleleSpecificAnnotationData<>(vc.getAlleles(),null);
            likelihoods.bestAllelesBreakingTies(sample).stream()
                    .filter(ba -> ba.isInformative())
                    .forEach(ba -> updateTable(ba.allele, ba.second_best_allele, ba.evidence, ref, allAlts, sampleTable, ba.likelihood, ba.secondBestLikelihood ));
            if (passesMinimumThreshold(sampleTable, minCount)) {
                combineAttributeMap(sampleTable, perAlleleValues);
            }
        }
    }


    private void updateTable(final Allele bestAllele, final Allele secondBestAllele, final GATKRead read, final Allele ref, final List<Allele> allAlts, final ReducibleAnnotationData<List<Double>> perAlleleValues,
                             final double lik_best, final double lik_second_best) {

        final boolean matchesRef = bestAllele.equals(ref, true);
        final boolean matchesAnyAlt = allAlts.contains(bestAllele);

        final boolean secondMatchesRef = secondBestAllele.equals(ref, true);
        final boolean secondMatchesAnyAlt = allAlts.contains(secondBestAllele);


        //can happen if a read's most likely allele has been removed when --max_alternate_alleles is exceeded
        if (!( matchesRef || matchesAnyAlt )) {
            return;
        }

        double nf = Math.pow(10, lik_best) + Math.pow(10, lik_second_best);
        double allocation_best = Math.pow(10, lik_best)/nf;
        double allocation_second_best = Math.pow(10, lik_second_best)/nf;

        if ((Double.isNaN(allocation_second_best))) {
            if ((Double.isNaN(allocation_best))) {
                return;
            }
            else {
                allocation_best = 1;
                allocation_second_best = 0 ;
            }

        }

        ArrayList<Allele> alleleList = new ArrayList<>();
        alleleList.add(0, bestAllele);
        alleleList.add(1, secondBestAllele);

        ArrayList<Double> allocationList = new ArrayList<>(2);
        allocationList.add(0, allocation_best);
        allocationList.add(1, allocation_second_best);


        for (int i = 0 ; i < alleleList.size(); i++ ) {
            Allele allele = alleleList.get(i);
            Double allocation = allocationList.get(i);
            List<Double> alleleStrandCounts;
            if (perAlleleValues.hasAttribute(allele) && perAlleleValues.getAttribute(allele) != null) {
                alleleStrandCounts = perAlleleValues.getAttribute(allele);
            } else {
                alleleStrandCounts = new ArrayList<>();
                alleleStrandCounts.add(0, (double) 0);
                alleleStrandCounts.add(1, (double) 0);
            }

            boolean isForward = !read.isReverseStrand();
            if (isForward) {
                alleleStrandCounts.set(FORWARD, alleleStrandCounts.get(FORWARD) + allocation);
            } else {
                alleleStrandCounts.set(REVERSE, alleleStrandCounts.get(REVERSE) + allocation);
            }
            perAlleleValues.putAttribute(allele, alleleStrandCounts);
        }
    }


    protected boolean passesMinimumThreshold(final ReducibleAnnotationData<List<Double>> sampleTable, final double minCount) {
        final double readCount = sampleTable.getAttributeMap().values().stream()
                .filter(alleleValues -> alleleValues != null)
                .mapToDouble(alleleValues -> alleleValues.get(FORWARD) + alleleValues.get(REVERSE))
                .sum();
        return readCount > minCount;
    }


    protected void combineAttributeMap(final ReducibleAnnotationData<List<Double>> toAdd, final ReducibleAnnotationData<List<Double>> combined) {
        for (final Allele a : combined.getAlleles()) {
            if (toAdd.hasAttribute(a) && toAdd.getAttribute(a) != null) {
                if (combined.getAttribute(a) != null) {
                    combined.getAttribute(a).set(FORWARD, combined.getAttribute(a).get(FORWARD) + toAdd.getAttribute(a).get(FORWARD));
                    combined.getAttribute(a).set(REVERSE, combined.getAttribute(a).get(REVERSE) + toAdd.getAttribute(a).get(REVERSE));
                }
                else {
                    List<Double> alleleData = new ArrayList<>();
                    alleleData.add(FORWARD, toAdd.getAttribute(a).get(FORWARD));
                    alleleData.add(REVERSE, toAdd.getAttribute(a).get(REVERSE));
                    combined.putAttribute(a,alleleData);
                }
            }
        }
    }

    @Override
    //Allele-specific annotations cannot be called from walkers other than HaplotypeCaller
    protected Map<String, Object> calculateAnnotationFromGTfield(final GenotypesContext genotypes){
        return Collections.emptyMap();
    }

    @Override
    //Allele-specific annotations cannot be called from walkers other than HaplotypeCaller
    protected Map<String, Object> calculateAnnotationFromStratifiedContexts(final Map<String, List<PileupElement>> stratifiedContexts,
                                                                            final VariantContext vc){
        return new HashMap<>();
    }

    public static String rawValueAsString(int[][] table) {
        return table[0][0]+","+table[0][1]+ PRINT_DELIM +table[1][0]+","+table[1][1];
    }

    /**
     * Method which combines the per allele contingency tables from the underlying variant contexts by totaling
     * supported values for both forward and reverse data and outputting it as a new contingency table.
     *
     * @param vcAlleles
     * @param annotationList
     * @return
     */
    @Override
    @SuppressWarnings({"unchecked", "rawtypes"})//FIXME generics here blow up
    public Map<String, Object> combineRawData(final List<Allele> vcAlleles, final List<ReducibleAnnotationData<?>>  annotationList) {
        //VC already contains merged alleles from ReferenceConfidenceVariantContextMerger
        ReducibleAnnotationData combinedData = new AlleleSpecificAnnotationData(vcAlleles, null);

        for (final ReducibleAnnotationData currentValue : annotationList) {
            parseRawDataString(currentValue);
            combineAttributeMap(currentValue, combinedData);
        }
        final String annotationString = makeRawAnnotationString(vcAlleles, combinedData.getAttributeMap());
        return Collections.singletonMap(getPrimaryRawKey(), annotationString);
    }

    /**
     * Parses the raw data stings of combined contingency matrix data and calls client methods calculateReducedData(myData)
     * implementation to generate double digest of provided allele information which is stored in '|' delineated lists.
     *
     * @param vc -- contains the final set of alleles, possibly subset by GenotypeGVCFs
     * @param originalVC -- used to get all the alleles for all gVCFs
     * @return
     */
    @Override
    public  Map<String, Object> finalizeRawData(final VariantContext vc, final VariantContext originalVC) {
        if (!vc.hasAttribute(getPrimaryRawKey())) {
            return new HashMap<>();
        }
        String rawRankSumData = vc.getAttributeAsString(getPrimaryRawKey(),null);
        if (rawRankSumData == null) {
            return new HashMap<>();
        }
        AlleleSpecificAnnotationData<List<Double>> myData = new AlleleSpecificAnnotationData<>(originalVC.getAlleles(), rawRankSumData);
        parseRawDataString(myData);

        Map<Allele, Double> perAltRankSumResults = calculateReducedData(myData);

        String annotationString = makeReducedAnnotationString(vc, perAltRankSumResults);
        return Collections.singletonMap(getKeyNames().get(0), annotationString);
    }

    protected String makeReducedAnnotationString(VariantContext vc, Map<Allele,Double> perAltsStrandCounts) {
        String annotationString = "";
        for (Allele a : vc.getAlternateAlleles()) {
            if (!annotationString.isEmpty()) {
                annotationString += REDUCED_DELIM;
            }
            if (!perAltsStrandCounts.containsKey(a)) {
                logger.warn("ERROR: VC allele not found in annotation alleles -- maybe there was trimming?");
            } else {
                annotationString += String.format("%.3f", perAltsStrandCounts.get(a));
            }
        }
        return annotationString;
    }


    protected void parseRawDataString(ReducibleAnnotationData<List<Double>> myData) {
        String rawDataString = myData.getRawData();
        if (rawDataString.startsWith("[")) {
            rawDataString = rawDataString.substring(1,rawDataString.length()-1);
        }
        String[] rawDataPerAllele;
        String[] rawListEntriesAsStringVector;
        Map<Allele, List<Double>> perAlleleValues = new HashMap<>();
        //Initialize maps
        for (Allele current : myData.getAlleles()) {
            perAlleleValues.put(current, new LinkedList<Double>());
        }
        //rawDataPerAllele is the list of values for each allele (each of variable length)
        rawDataPerAllele = rawDataString.split(SPLIT_DELIM);
        for (int i=0; i<rawDataPerAllele.length; i++) {
            String alleleData = rawDataPerAllele[i];
            if (!alleleData.isEmpty()) {
                List<Double> alleleList = perAlleleValues.get(myData.getAlleles().get(i));
                rawListEntriesAsStringVector = alleleData.split(",");
                //Read counts will only ever be integers
                for (String s : rawListEntriesAsStringVector) {
                    if (!s.isEmpty()) {
                        alleleList.add(Double.parseDouble(s.trim()));
                    }
                }
            }
        }
        myData.setAttributeMap(perAlleleValues);
    }

    /**
     * Method which determines how the Strand Bias read direction allele data must be combined into a final annotation
     * Must be overridden by client methods.
     *
     * @param combinedData
     * @return
     */
    protected abstract Map<Allele,Double> calculateReducedData(final AlleleSpecificAnnotationData<List<Double>> combinedData );

}
