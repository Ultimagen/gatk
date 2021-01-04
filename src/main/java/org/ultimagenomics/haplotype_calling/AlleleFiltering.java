package org.ultimagenomics.haplotype_calling;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.InverseAllele;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.ultimagenomics.flow_based_read.read.FlowBasedHaplotype;

import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public abstract class AlleleFiltering {

    protected static final Logger logger = LogManager.getLogger(AlleleFiltering.class);
    private AssemblyBasedCallerArgumentCollection hcArgs;
    private PrintStream assemblyDebugOutStream;

    public AlleleFiltering(AssemblyBasedCallerArgumentCollection _hcargs, PrintStream _assemblyDebugOutStream){
        hcArgs = _hcargs;
        assemblyDebugOutStream = _assemblyDebugOutStream;
    }

    public AlleleLikelihoods<GATKRead, Haplotype> filterAlleles(AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods, int activeWindowStart){
        AlleleLikelihoods<GATKRead, Haplotype> subsettedReadLikelihoodsFinal;

        logger.debug("SHA:: filter alleles - start");
        AlleleLikelihoods<GATKRead, Haplotype> subsettedReadLikelihoodsBoth = subsetHaplotypesByAlleles(readLikelihoods, hcArgs, activeWindowStart);
        logger.debug("SHA:: filter alleles - end");

        subsettedReadLikelihoodsFinal = subsettedReadLikelihoodsBoth;
        if (assemblyDebugOutStream != null) {
            assemblyDebugOutStream.println("\nThere were " + subsettedReadLikelihoodsFinal.alleles().size() + " haplotypes found after subsetting by alleles. Here they are:");
            subsettedReadLikelihoodsFinal.alleles().stream().map(Haplotype::toString).sorted().forEach(assemblyDebugOutStream::println);
        }

        return subsettedReadLikelihoodsFinal;
    }



    static private Set<LocationAndAllele> getAlleles(final Haplotype haplotype){
        Collection<VariantContext> vcs = haplotype.getEventMap().getVariantContexts();
        Set<LocationAndAllele> allEvents = new HashSet<>();
        for (VariantContext vc: vcs) {
            allEvents.addAll(vc.getAlleles().stream().map( al -> new LocationAndAllele(vc.getStart(), al, vc.getReference())).collect(Collectors.toList()));
        }
        return allEvents;
    }

    private AlleleLikelihoods<GATKRead, Haplotype> subsetHaplotypesByAlleles(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                                             AssemblyBasedCallerArgumentCollection hcargs,
                                                                             int activeWindowStart){

        boolean removedHaplotype = true;
        AlleleLikelihoods<GATKRead, Haplotype> currentReadLikelihoods = readLikelihoods;
        final Map<Haplotype, Collection<LocationAndAllele>> haplotypeAlleleMap  = new CollectionUtil.DefaultingMap<>((k) -> new ArrayList<>(), true);
        readLikelihoods.alleles().forEach(h -> getAlleles(h).stream().filter(al -> !al.isReference()).forEach(jh -> haplotypeAlleleMap.get(h).add(jh)));
        OccurrenceMatrix<Haplotype, LocationAndAllele> occm = new OccurrenceMatrix<>(haplotypeAlleleMap);
        List<Pair<LocationAndAllele, LocationAndAllele>> nonCoOcurringAlleles = occm.nonCoOcurringColumns();
        nonCoOcurringAlleles = filterByDistance(nonCoOcurringAlleles, 1, 20);
        nonCoOcurringAlleles = filterSameUpToHmerPairs(nonCoOcurringAlleles, findReferenceHaplotype(readLikelihoods.alleles()), activeWindowStart);

        List<Set<LocationAndAllele>> independentAlleles = occm.getIndependentSets(nonCoOcurringAlleles);
        List<LocationAndAllele> allRemovedAlleles = new ArrayList<>();
        Set<Haplotype> haplotypesToRemove = new HashSet<>();


        for ( Set<LocationAndAllele> alleleSet : independentAlleles) {
            Set<Haplotype> enabledHaplotypes = new HashSet<>();

            for (Haplotype h : currentReadLikelihoods.alleles()) enabledHaplotypes.add(h);

            while (removedHaplotype) {
                // build map from contig to haplotype
                final Map<LocationAndAllele, List<Haplotype>> alleleHaplotypeMap = new CollectionUtil.DefaultingMap<>((k) -> new ArrayList<>(), true);

                currentReadLikelihoods.alleles().stream().filter(h->enabledHaplotypes.contains(h))
                        .forEach(h -> getAlleles(h).stream().filter(al -> alleleSet.contains(al)).filter(al -> !al.isReference()).forEach(
                                jh -> alleleHaplotypeMap.get(jh).add(h))
                        );


                logger.debug("AHM::printout start");
                for (LocationAndAllele al : alleleHaplotypeMap.keySet()) {
                    logger.debug("AHM::allele block ---> ");
                    for (Allele h : alleleHaplotypeMap.get(al)) {
                        logger.debug(String.format("AHM:: (%d) %s: %s", al.getLoc(), al.getAllele().getBaseString(), h.getBaseString()));
                    }
                    logger.debug("AHM::allele block ---< ");

                }
                logger.debug("AHM::printout end");

                final List<Haplotype> eventualAlleles = new ArrayList<>(currentReadLikelihoods.alleles());
                if (eventualAlleles.stream().noneMatch(Allele::isReference)) {
                    throw new IllegalStateException("Reference haplotype must always remain!");
                }
                // repeat until no change.
                removedHaplotype = false;


                // find all the haplotypes have them and remove them from a consideration to be removed.
                final List<LocationAndAllele> allHapAlleles = alleleHaplotypeMap.keySet().stream().filter(c -> alleleHaplotypeMap.get(c).containsAll(eventualAlleles)).collect(Collectors.toList());
                allHapAlleles.forEach(alleleHaplotypeMap::remove);

                logger.debug("----- AHA start ----");
                for (LocationAndAllele al : allHapAlleles) {
                    logger.debug(String.format("AHA:: %s", (al).toString()));
                }
                logger.debug("----- AHA end -----");


                //find alleles that no haplotypes have them and remove them from a consideration to be removed (nobody cares).
                final List<LocationAndAllele> noHapAlleles = alleleHaplotypeMap.keySet().stream().filter(c -> alleleHaplotypeMap.get(c).isEmpty()).collect(Collectors.toList());
                noHapAlleles.forEach(alleleHaplotypeMap::remove);
                logger.debug("----- NHA start ----");
                for (LocationAndAllele al : noHapAlleles) {
                    logger.debug(String.format("NHA:: %s", (al.toString())));
                }
                logger.debug("----- NHA end -----");

                final List<LocationAndAllele> allAlleles = new ArrayList<>(alleleHaplotypeMap.keySet());

                final AlleleLikelihoods<GATKRead, Haplotype> finalCurrentReadLikelihoods = currentReadLikelihoods;
                logger.debug("GAL::start");
                final List<AlleleLikelihoods<GATKRead, Allele>> alleleLikelihoods =
                        allAlleles.stream().map(c -> getAlleleLikelihoodMatrix(finalCurrentReadLikelihoods, c, haplotypeAlleleMap)).collect(Collectors.toList());

                final List<Integer> collectedRPLs = IntStream.range(0, allAlleles.size()).mapToObj(i -> getAlleleLikelihood(alleleLikelihoods.get(i), allAlleles.get(i))).collect(Collectors.toList());
                final List<Double> collectedSORs = IntStream.range(0, allAlleles.size()).mapToObj(i -> getAlleleSOR(alleleLikelihoods.get(i), allAlleles.get(i))).collect(Collectors.toList());
                logger.debug("GAL::end");

                final LocationAndAllele badAllele;
                badAllele = identifyBadAllele(collectedRPLs, collectedSORs, allAlleles, hcargs.prefilterQualThreshold, hcargs.prefilterSorThreshold);


                if (badAllele != null) {
                    logger.debug(String.format("SHA:: Removing %s", badAllele.toString()));
                    allRemovedAlleles.add(badAllele);
                    removedHaplotype = true;
                    final ArrayList<Haplotype> haplotypesWithAllele = new ArrayList<>(alleleHaplotypeMap.get(badAllele));
                    haplotypesWithAllele.removeIf(Allele::isReference);
                    for (Haplotype h: haplotypesWithAllele){
                        haplotypesToRemove.add(h);
                        enabledHaplotypes.remove(h);
                    }
                }
            }
        }

        Set<Haplotype> eventualAlleles = new HashSet<>();
        currentReadLikelihoods.alleles().stream().filter(al -> !haplotypesToRemove.contains(al)).forEach(al -> eventualAlleles.add(al));
        currentReadLikelihoods = currentReadLikelihoods.subsetToAlleles(eventualAlleles);

        return currentReadLikelihoods;
    }


    private LocationAndAllele identifyBadAllele(final List<Integer> collectedRPLs,final List<Double> collectedSORs, List<LocationAndAllele> alleles, double qual_threshold, double sor_threshold) {
        final Integer worstRPL = collectedRPLs.stream().max(Integer::compareTo).orElse(Integer.MIN_VALUE);
        final Double worstSOR = collectedSORs.stream().max(Double::compareTo).orElse(0.0);
        final double THRESHOLD = -1 * qual_threshold;
        final double SOR_THRESHOLD = sor_threshold;

        int argmin=-1;
        if (worstRPL > THRESHOLD ) {
            argmin = collectedRPLs.indexOf(worstRPL);
            logger.debug(String.format("SHA:: WorstRPL: %d", worstRPL));
        }
        else if (worstSOR > SOR_THRESHOLD ) {
            argmin = collectedSORs.indexOf(worstSOR);
            logger.debug(String.format("SHA:: WorstSOR: %f", worstSOR));

        }
        if (argmin >= 0 ) {
            logger.debug(String.format("SHA:: %d", argmin));
            return alleles.get(argmin);
        } else
            return null;
    }


    private enum StrandsToUse {
        FORWARD(true),
        REVERSE(false);

        private final boolean useForwardStrand;

        StrandsToUse(final boolean useForwardStrand) {
            this.useForwardStrand = useForwardStrand;
        }

        boolean useRead(GATKRead read) {
            return read.isReverseStrand() ^ useForwardStrand;
        }
    }


    private AlleleLikelihoods<GATKRead, Allele> getAlleleLikelihoodMatrix(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                                                     final LocationAndAllele allele,
                                                                                     final Map<Haplotype, Collection<LocationAndAllele>> haplotypeAlleleMap
                                                                                     ){
        Map<Allele,List<Haplotype>> alleleHaplotypeMap = new CollectionUtil.DefaultingMap<>((k) -> new ArrayList<>(), true);

        final Allele notAllele= InverseAllele.of(allele.getAllele());
        List<Allele> tmp = new ArrayList<>(1);
        tmp.add(allele);


        readLikelihoods.alleles().stream().filter(h->haplotypeAlleleMap.get(h).contains(allele)).forEach(alleleHaplotypeMap.get(allele)::add);
        readLikelihoods.alleles().stream().filter(h -> !haplotypeAlleleMap.get(h).contains(allele)).forEach(alleleHaplotypeMap.get(notAllele)::add);

        final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods = readLikelihoods.marginalize(alleleHaplotypeMap);
        logger.debug(String.format("GALM: %s %d %d", allele.toString(), alleleHaplotypeMap.get(allele).size(), alleleHaplotypeMap.get(notAllele).size()));
        return alleleLikelihoods;
    }


    abstract int getAlleleLikelihood(final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods, Allele allele);
    abstract double getAlleleSOR(final AlleleLikelihoods<GATKRead, Allele> alleleLikelihoods, Allele allele);


    static public List<Pair<LocationAndAllele, LocationAndAllele>> filterByDistance(List<Pair<LocationAndAllele, LocationAndAllele>> allelePairs,
                                   int min_dist, int max_dist) {
        logger.debug(String.format("FBD: input %d pairs ", allelePairs.size()));
        allelePairs.removeIf(v -> Math.abs(v.getLeft().getLoc() - v.getRight().getLoc())>max_dist);
        allelePairs.removeIf(v -> Math.abs(v.getLeft().getLoc() - v.getRight().getLoc())<min_dist);
        logger.debug(String.format("FBD: output %d pairs ", allelePairs.size()));

        return allelePairs;
    }


    public List<Pair<LocationAndAllele, LocationAndAllele>> filterSameUpToHmerPairs(List<Pair<LocationAndAllele,
            LocationAndAllele>> allelePairs, Haplotype refHaplotype, int activeWindowStart) {
        logger.debug(String.format("FSUHP: input %d pairs ", allelePairs.size()));

        List<Pair<LocationAndAllele, LocationAndAllele>> result = new ArrayList<>();

        for (Pair<LocationAndAllele, LocationAndAllele> allelePair: allelePairs) {

            int commonPrefixLengthLeft = Math.min(allelePair.getLeft().getAllele().length(), allelePair.getLeft().getRefAllele().length());
            int commonPrefixLengthRight = Math.min(allelePair.getRight().getAllele().length(), allelePair.getRight().getRefAllele().length());

            Pair<Haplotype, Haplotype> modifiedHaplotypes = new ImmutablePair<>(
                    refHaplotype.insertAllele(
                            allelePair.getLeft().getRefAllele(),
                            allelePair.getLeft().getAllele(),
                            allelePair.getLeft().getLoc()-activeWindowStart,
                            -1,
                            commonPrefixLengthLeft),
                    refHaplotype.insertAllele(
                            allelePair.getRight().getRefAllele(),
                            allelePair.getRight().getAllele(),
                            allelePair.getRight().getLoc() - activeWindowStart,
                            -1,
                            commonPrefixLengthRight));
            logger.debug("FSUHP: Insertion:");
            logger.debug(String.format("FSUHP: refAllele: (%d) %s", refHaplotype.getStartPosition(), refHaplotype.toString()));
            logger.debug(String.format("FSUHP: Replace: %s", allelePair.getRight().toString()));
            logger.debug(String.format("FSUHP: result: %s", modifiedHaplotypes.getRight().toString()));

            Pair<FlowBasedHaplotype, FlowBasedHaplotype> fbh = new ImmutablePair<>(haplotype2FlowHaplotype(modifiedHaplotypes.getLeft()),
                                                                                   haplotype2FlowHaplotype(modifiedHaplotypes.getRight()));

            if (fbh.getLeft().equalUpToHmerChange(fbh.getRight()))
            {
                result.add(allelePair);
            }
        }

        logger.debug(String.format("FSUHP: output %d pairs ", result.size()));

        return result;
    }


    private Haplotype findReferenceHaplotype( List<Haplotype> haplotypeList) {
        for (Haplotype h: haplotypeList ) {
            if (h.isReference()) {
                return h;
            }
        }
        return null;
    }


    private FlowBasedHaplotype haplotype2FlowHaplotype(Haplotype hap) {
        FlowBasedHaplotype flowBasedHaplotype = new FlowBasedHaplotype(hap, "TACG");
        return flowBasedHaplotype;
    }

}
