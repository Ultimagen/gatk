package org.ultimagenomics.haplotype_calling;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.InverseAllele;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

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

    public AlleleLikelihoods<GATKRead, Haplotype> filterAlleles(AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods){
        AlleleLikelihoods<GATKRead, Haplotype> subsettedReadLikelihoodsFinal;

        logger.debug("SHA:: filter alleles - start");
        AlleleLikelihoods<GATKRead, Haplotype> subsettedReadLikelihoodsBoth = subsetHaplotypesByAlleles(readLikelihoods, hcArgs);
        logger.debug("SHA:: filter alleles - end");

        subsettedReadLikelihoodsFinal = subsettedReadLikelihoodsBoth;
        if (assemblyDebugOutStream != null) {
            assemblyDebugOutStream.println("\nThere were " + subsettedReadLikelihoodsFinal.alleles().size() + " haplotypes found after subsetting by alleles. Here they are:");
            subsettedReadLikelihoodsFinal.alleles().stream().map(Haplotype::toString).sorted().forEach(assemblyDebugOutStream::println);
        }

        return subsettedReadLikelihoodsFinal;
    }

    private AlleleLikelihoods<GATKRead, Haplotype> subsetHaplotypesByAlleles(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                                             AssemblyBasedCallerArgumentCollection hcargs) {
        return subsetHaplotypesByAlleles(readLikelihoods, hcargs, true);
    }


    static private Set<LocationAndAllele> getAlleles(final Haplotype haplotype){
        Collection<VariantContext> vcs = haplotype.getEventMap().getVariantContexts();
        Set<LocationAndAllele> allEvents = new HashSet<>();
        for (VariantContext vc: vcs) {
            allEvents.addAll(vc.getAlleles().stream().map( al -> new LocationAndAllele(vc.getStart(), al)).collect(Collectors.toList()));
        }
        return allEvents;
    }

    private AlleleLikelihoods<GATKRead, Haplotype> subsetHaplotypesByAlleles(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                                             AssemblyBasedCallerArgumentCollection hcargs, boolean keepRef ){

        boolean removedHaplotype = true;
        AlleleLikelihoods<GATKRead, Haplotype> currentReadLikelihoods = readLikelihoods;
        final Map<Haplotype, List<LocationAndAllele>> haplotypeAlleleMap  = new CollectionUtil.DefaultingMap<>((k) -> new ArrayList<>(), true);
        readLikelihoods.alleles().forEach(h -> getAlleles(h).forEach(jh -> haplotypeAlleleMap.get(h).add(jh)));

        while (removedHaplotype) {
            // build map from contig to haplotype
            final Map<LocationAndAllele, List<Haplotype>> alleleHaplotypeMap = new CollectionUtil.DefaultingMap<>((k) -> new ArrayList<>(), true);

            currentReadLikelihoods.alleles()
                    .forEach(h -> getAlleles(h).forEach(
                            jh -> alleleHaplotypeMap.get(jh).add(h))
                    );


            logger.debug("AHM::printout start");
            for (LocationAndAllele al: alleleHaplotypeMap.keySet()) {
                logger.debug("AHM::allele block ---> ");
                for (Allele h: alleleHaplotypeMap.get(al)){
                    logger.debug(String.format("AHM:: (%d) %s: %s", al.getLoc(), al.getAllele().getBaseString(),h.getBaseString()));
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
            for (LocationAndAllele al: allHapAlleles) {
                logger.debug(String.format("AHA:: %s", (al).toString()));
            }
            logger.debug("----- AHA end -----");


            //find alleles that no haplotypes have them and remove them from a consideration to be removed (nobody cares).
            final List<LocationAndAllele> noHapAlleles= alleleHaplotypeMap.keySet().stream().filter(c -> alleleHaplotypeMap.get(c).isEmpty()).collect(Collectors.toList());
            noHapAlleles.forEach(alleleHaplotypeMap::remove);
            logger.debug("----- NHA start ----");
            for (LocationAndAllele al: noHapAlleles) {
                logger.debug(String.format("NHA:: %s", (al.toString())));
            }
            logger.debug("----- NHA end -----");

            final List<LocationAndAllele> allAlleles = new ArrayList<>(alleleHaplotypeMap.keySet());

            final AlleleLikelihoods<GATKRead, Haplotype> finalCurrentReadLikelihoods = currentReadLikelihoods;
            logger.debug("GAL::start");
            final List<AlleleLikelihoods<GATKRead, Allele>> alleleLikelihoods =
                    allAlleles.stream().map(c -> getAlleleLikelihoodMatrix(finalCurrentReadLikelihoods, c, haplotypeAlleleMap)).collect(Collectors.toList());

            final List<Integer> collectedRPLs = IntStream.range(0, allAlleles.size()).mapToObj(i -> getAlleleLikelihood(alleleLikelihoods.get(i), allAlleles.get(i))).collect(Collectors.toList());
            final List<Double> collectedSORs = IntStream.range(0, allAlleles.size()).mapToObj( i -> getAlleleSOR(alleleLikelihoods.get(i), allAlleles.get(i))).collect(Collectors.toList());
            logger.debug("GAL::end");

            final LocationAndAllele badAllele;
            badAllele = identifyBadAllele(collectedRPLs, collectedSORs, allAlleles, hcargs.prefilterQualThreshold, hcargs.prefilterSorThreshold);


            if (badAllele != null){
                logger.debug(String.format("SHA:: Removing %s", badAllele.toString()));
                final ArrayList<Haplotype> haplotypesWithAllele = new ArrayList<>(alleleHaplotypeMap.get(badAllele));
                haplotypesWithAllele.removeIf(Allele::isReference);

                removedHaplotype = eventualAlleles.removeAll(haplotypesWithAllele);
                if (eventualAlleles.stream().noneMatch(Allele::isReference)) {
                    throw new IllegalStateException("Reference haplotype must always remain!");
                }
                //subset to remaining haplotypes:
                currentReadLikelihoods = currentReadLikelihoods.subsetToAlleles(eventualAlleles);
            }
        }

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
                                                                                     final Map<Haplotype, List<LocationAndAllele> > haplotypeAlleleMap
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

}
