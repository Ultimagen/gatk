package org.ultimagenomics.haplotype_calling;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.variant.variantcontext.Allele;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.broadinstitute.hellbender.tools.walkers.annotator.StrandOddsRatio;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingData;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingLikelihoods;
import org.broadinstitute.hellbender.tools.walkers.genotyper.IndependentSampleGenotypesModel;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.InverseAllele;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.JoinedContigs;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ContigFiltering {
    private static final Logger logger = LogManager.getLogger(ContigFiltering.class);
    private HaplotypeCallerArgumentCollection hcArgs;
    public ContigFiltering(HaplotypeCallerArgumentCollection _hcargs){
        hcArgs = _hcargs;
    }
    public AlleleLikelihoods<GATKRead, Haplotype> filterContigs(AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods){
        AlleleLikelihoods<GATKRead, Haplotype> subsettedReadLikelihoodsFinal;

        logger.debug("SHC:: filter contigs - start");
        AlleleLikelihoods<GATKRead, Haplotype> subsettedReadLikelihoodsBoth = subsetHaplotypesByContigs(readLikelihoods, hcArgs);
        logger.debug("SHC:: filter contigs - end");
        logger.debug("SHC:: filter ref contigs - start");
        subsettedReadLikelihoodsBoth = subsetHaplotypesByRefContigs(subsettedReadLikelihoodsBoth, hcArgs);
                logger.debug("SHC:: filter ref contigs - end");

        subsettedReadLikelihoodsFinal = subsettedReadLikelihoodsBoth;
        if (assemblyDebugOutStream != null) {
            assemblyDebugOutStream.println("\nThere were " + subsettedReadLikelihoodsFinal.alleles().size() + " haplotypes found after subsetting by contigs. Here they are:");
            subsettedReadLikelihoodsFinal.alleles().stream().map(Haplotype::toString).sorted().forEach(assemblyDebugOutStream::println);
        }

        return subsettedReadLikelihoodsFinal;
    }

    private AlleleLikelihoods<GATKRead, Haplotype> subsetHaplotypesByContigs(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                                             HaplotypeCallerArgumentCollection hcargs) {
        return subsetHaplotypesByContigs(readLikelihoods, hcargs, true);
    }

    private AlleleLikelihoods<GATKRead, Haplotype> subsetHaplotypesByRefContigs(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                                                HaplotypeCallerArgumentCollection hcargs) {
        return subsetHaplotypesByContigs(readLikelihoods, hcargs, false);
    }

    static private Set<Allele> getJoinedContigs(final Haplotype haplotype){
        Set<Allele> joinedContigs = new HashSet<>();
        //noinspection ResultOfMethodCallIgnored
        haplotype.contigs.stream().reduce((a, b) -> {
            joinedContigs.add(new JoinedContigs(a, b));
            return b;
        });
        return joinedContigs;
    }

    private AlleleLikelihoods<GATKRead, Haplotype> subsetHaplotypesByContigs(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods,
                                                                             HaplotypeCallerArgumentCollection hcargs, boolean keepRef ){

        boolean removedHaplotype = true;
        AlleleLikelihoods<GATKRead, Haplotype> currentReadLikelihoods = readLikelihoods;
        while (removedHaplotype) {
            // build map from contig to haplotype
            final Map<Allele, List<Haplotype>> contigHaplotypeMap = new CollectionUtil.DefaultingMap<>((k) -> new ArrayList<>(), true);

            currentReadLikelihoods.alleles()
                    .forEach(h -> getJoinedContigs(h).forEach(
                            jh -> contigHaplotypeMap.get(jh).add(h))
                    );


            logger.debug("CHM::printout start");
            for (Allele jc: contigHaplotypeMap.keySet()) {
                logger.debug("CHM::contig block ---> ");
                for (Allele h: contigHaplotypeMap.get(jc)){
                    logger.debug(String.format("CHM:: %s->%s: %s", ((JoinedContigs)jc).getAllele1().getBaseString(),
                            ((JoinedContigs)jc).getAllele2().getBaseString(), h.getBaseString()));
                }
                logger.debug("CHM::contig block ---< ");

            }
            logger.debug("CHM::printout end");

            final List<Haplotype> eventualAlleles = new ArrayList<>(currentReadLikelihoods.alleles());
            if (eventualAlleles.stream().noneMatch(Allele::isReference)) {
                throw new IllegalStateException("Reference haplotype must always remain!");
            }
            // repeat until no change.
            removedHaplotype = false;

            //find contigs that only have the reference haplotypes and remove them.
            if (keepRef) {
                final List<Allele> refOnlyContigs = contigHaplotypeMap.keySet().stream().filter(c -> contigHaplotypeMap.get(c).stream().anyMatch(Allele::isReference)).collect(Collectors.toList());
                refOnlyContigs.forEach(contigHaplotypeMap::remove);
                logger.debug("----- ROC start ");
                for (Allele al : refOnlyContigs) {
                    logger.debug(String.format("ROC:: %s", (al).toString()));
                }
                logger.debug("----- ROC end");
            }

            // find contigs that have all the haplotypes in them and remove them.
            final List<Allele> allHapContigs = contigHaplotypeMap.keySet().stream().filter(c -> contigHaplotypeMap.get(c).containsAll(eventualAlleles)).collect(Collectors.toList());
            allHapContigs.forEach(contigHaplotypeMap::remove);

            logger.debug("----- AHC start ----");
            for (Allele al: allHapContigs) {
                logger.debug(String.format("AHC:: %s", (al).toString()));
            }
            logger.debug("----- AHC end -----");


            //find contigs that have no haplotypes in them and remove them.
            final List<Allele> noHapContigs = contigHaplotypeMap.keySet().stream().filter(c -> contigHaplotypeMap.get(c).isEmpty()).collect(Collectors.toList());
            noHapContigs.forEach(contigHaplotypeMap::remove);
            logger.debug("----- NHC start ----");
            for (Allele al: noHapContigs) {
                logger.debug(String.format("NHC:: %s", (al.toString())));
            }
            logger.debug("----- NHC end -----");

            final List<Allele> allAlleles = new ArrayList<>(contigHaplotypeMap.keySet());

            final AlleleLikelihoods<GATKRead, Haplotype> finalCurrentReadLikelihoods = currentReadLikelihoods;
            logger.debug("GCL::start");
            final List<AlleleLikelihoods<GATKRead, Allele>> contigLikelihoods =
                    allAlleles.stream().map(c -> getContigLikelihoodMatrix(finalCurrentReadLikelihoods, c)).collect(Collectors.toList());

            final List<Integer> collectedRPLs = IntStream.range(0, allAlleles.size()).mapToObj(i -> getContigLikelihood(contigLikelihoods.get(i), allAlleles.get(i))).collect(Collectors.toList());
            final List<Double> collectedSORs = IntStream.range(0, allAlleles.size()).mapToObj( i -> getContigSOR(contigLikelihoods.get(i), allAlleles.get(i))).collect(Collectors.toList());
            logger.debug("GCL::end");

            final Allele badContig;
            if (keepRef) {
                badContig = identifyBadContig(collectedRPLs, collectedSORs, allAlleles, hcargs.prefilterQualThreshold, hcargs.prefilterSorThreshold);
            } else {
                badContig = identifyBadContig(collectedRPLs, collectedSORs, allAlleles, hcargs.refPrefilterQualThreshold, hcargs.prefilterSorThreshold);
            }

            if (badContig != null){
                logger.debug(String.format("SHC:: Removing %s (%d) -> (%d)", badContig.toString(),
                        ((JoinedContigs)badContig).getAllele1().hashCode(),
                        ((JoinedContigs)badContig).getAllele2().hashCode()));
                final ArrayList<Haplotype> haplotypesWithContig = new ArrayList<>(contigHaplotypeMap.get(badContig));
                haplotypesWithContig.removeIf(Allele::isReference);

                removedHaplotype = eventualAlleles.removeAll(haplotypesWithContig);
                if (eventualAlleles.stream().noneMatch(Allele::isReference)) {
                    throw new IllegalStateException("Reference haplotype must always remain!");
                }
                //subset to remaining haplotypes:
                currentReadLikelihoods = currentReadLikelihoods.subsetToAlleles(eventualAlleles);
            }
        }

        return currentReadLikelihoods;
    }

    private Allele identifyBadContig(final List<Integer> collectedRPLs,final List<Double> collectedSORs, List<Allele> alleles, double qual_threshold, double sor_threshold) {
        final Integer worstRPL = collectedRPLs.stream().max(Integer::compareTo).orElse(Integer.MIN_VALUE);
        final Double worstSOR = collectedSORs.stream().max(Double::compareTo).orElse(0.0);
        final double THRESHOLD = -1 * qual_threshold;
        final double SOR_THRESHOLD = sor_threshold;

        int argmin=-1;
        if (worstRPL > THRESHOLD ) {
            argmin = collectedRPLs.indexOf(worstRPL);
            logger.debug(String.format("SHC:: WorstRPL: %d", worstRPL));
        }
        else if (worstSOR > SOR_THRESHOLD ) {
            argmin = collectedSORs.indexOf(worstSOR);
            logger.debug(String.format("SHC:: WorstSOR: %f", worstSOR));

        }
        if (argmin >= 0 ) {
            logger.debug(String.format("SHC:: %d", argmin));
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

    private AlleleLikelihoods<GATKRead, Allele> getContigLikelihoodMatrix(final AlleleLikelihoods<GATKRead, Haplotype> readLikelihoods, final Allele contig){
        Map<Allele,List<Haplotype>> contigHaplotypeMap = new CollectionUtil.DefaultingMap<>((k) -> new ArrayList<>(), true);

        final Allele notContig = InverseAllele.of(contig);
        List<Allele> tmp = new ArrayList<>(1);
        tmp.add(contig);


        readLikelihoods.alleles().stream().filter(h -> getJoinedContigs(h).contains(contig)).forEach(contigHaplotypeMap.get(contig)::add);
        readLikelihoods.alleles().stream().filter(h -> !getJoinedContigs(h).contains(contig)).forEach(contigHaplotypeMap.get(notContig)::add);

        final AlleleLikelihoods<GATKRead, Allele> contigLikelihoods = readLikelihoods.marginalize(contigHaplotypeMap);
        logger.debug(String.format("GCLM: %s -> %s %d %d", ((JoinedContigs)contig).getAllele1().getBaseString(),
                ((JoinedContigs)contig).getAllele2().getBaseString(), contigHaplotypeMap.get(contig).size(), contigHaplotypeMap.get(notContig).size()));
        return contigLikelihoods;
    }


    private int getContigLikelihood(final AlleleLikelihoods<GATKRead, Allele> contigLikelihoods, Allele contig) {
        final Allele notContig = InverseAllele.of(contig);

        // iterate over contigs and see what their qual is.

        GenotypingData<Allele> genotypingData = new GenotypingData<>(genotypingEngine.getPloidyModel(), contigLikelihoods);

        IndependentSampleGenotypesModel genotypesModel = new IndependentSampleGenotypesModel();

        AlleleList<Allele> alleleList = new IndexedAlleleList<>(Arrays.asList(contig, notContig));

        final GenotypingLikelihoods<Allele> genotypingLikelihoods = genotypesModel.calculateLikelihoods(alleleList, genotypingData);

        final int[] asPL = genotypingLikelihoods.sampleLikelihoods(0).getAsPLs();
        final int retVal;
        retVal = Math.min(asPL[1], asPL[0]) - asPL[2]; // if this is "large", reject the contig.
        logger.debug(String.format("GCL:: %s->%s: %d %d %d", ((JoinedContigs)contig).getAllele1().getBaseString(),
                ((JoinedContigs)contig).getAllele2().getBaseString(), asPL[0], asPL[1], asPL[2]));
        return retVal;

    }

    private double getContigSOR(final AlleleLikelihoods<GATKRead, Allele> contigLikelihoods, Allele contig) {
        final Allele notContig = InverseAllele.of(contig);
        int [][] contingency_table = StrandOddsRatio.getContingencyTable(contigLikelihoods, notContig, Arrays.asList(contig), 1);
        double sor = StrandOddsRatio.calculateSOR(contingency_table);
        logger.debug(String.format("GCS:: %s->%s: %f (%d %d %d %d)", ((JoinedContigs)contig).getAllele1().getBaseString(),
                ((JoinedContigs)contig).getAllele2().getBaseString(), sor, contingency_table[0][0], contingency_table[0][1], contingency_table[1][0], contingency_table[1][1]));
        return sor;

    }


}
