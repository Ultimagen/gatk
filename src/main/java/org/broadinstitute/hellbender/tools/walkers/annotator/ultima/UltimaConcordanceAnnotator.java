package org.broadinstitute.hellbender.tools.walkers.annotator.ultima;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.GenotypeAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardMutectAnnotation;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.ultimagen.flowBasedRead.read.FlowBasedRead;

import java.util.*;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="UltimaConcordanceAnnotator")
public class UltimaConcordanceAnnotator extends GenotypeAnnotation implements StandardMutectAnnotation {
    private final static Logger logger = LogManager.getLogger(UltimaConcordanceAnnotator.class);

    // additional constants
    protected static final String   C_INSERT = "ins";
    protected static final String   C_DELETE = "del";
    protected static final String   C_CSS_NA = "NA";
    protected static final String   C_CSS_CS = "cycle-skip";
    protected static final String   C_CSS_PCS = "possible-cycle-skip";
    protected static final String   C_CSS_NS = "non-skip";

    protected static final int      MOTIF_SIZE = 5;
    protected static final int      GC_CONTENT_SIZE = 10;
    protected static final String   FLOW_ORDER_DEFAULT = "TACG";

    static class LocalAttributes {
        ReferenceContext ref;
        GenotypeBuilder gb;
        String      flowOrder;

        String      indel;
        int         indelLength;
        int         hmerIndelLength;
        String      leftMotif;
        String      rightMotif;
    }

    @Override
    public void annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        Utils.nonNull(ref);
        Utils.nonNull(vc);
        Utils.nonNull(gb);

        // some annotators share results
        LocalAttributes         la = new LocalAttributes();
        la.ref = ref;
        la.gb = gb;

        // establish flow order
        la.flowOrder = establishFlowOrder(likelihoods);

        // call annotatotrs
        indelClassify(vc, la);
        isHmerIndel(vc, la);
        getMotif(vc, la);
        gcContent(vc, la);
        cycleSkip(vc, la);
    }

    private String establishFlowOrder(AlleleLikelihoods<GATKRead, Allele> likelihoods) {

        // extract from a read
        if ( likelihoods != null ) {
            List<GATKRead>  reads = likelihoods.sampleEvidence(0);
            if ( reads.size() > 0 ) {
                GATKRead        read = reads.get(0);
                if ( read instanceof FlowBasedRead ) {
                    return ((FlowBasedRead)read).getFlowOrder();
                }
            }
        }

        // if here, return default
        return FLOW_ORDER_DEFAULT;
    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(
                GATKVCFConstants.ULTIMA_INDEL_CLASSIFY, GATKVCFConstants.ULTIMA_INDEL_LENGTH,
                GATKVCFConstants.ULTIMA_HMER_INDEL_LENGTH, GATKVCFConstants.ULTIMA_HMER_INDEL_NUC,
                GATKVCFConstants.ULTIMA_LEFT_MOTIF, GATKVCFConstants.ULTIMA_RIGHT_MOTIF,
                GATKVCFConstants.ULTIMA_GC_CONTENT,
                GATKVCFConstants.ULTIMA_CYCLESKIP_STATUS);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return getKeyNames().stream().map(s -> GATKVCFHeaderLines.getFormatLine(s)).collect(Collectors.toList());
    }

    // "indel_classify" and "indel_length"
    private void indelClassify(final VariantContext vc, final LocalAttributes la) {

        if ( vc.isIndel() ) {

            /*
            if not x['indel']:
                return None
            elif len(x['ref']) < max([len(y) for y in x['alleles']]):
                return 'ins'
            return 'del'
             */
            final int maxAlleleLength = vc.getAlleles().stream()
                    .filter(allele -> !allele.isReference())
                    .map(allele -> allele.length())
                    .max(Integer::compare).get();
            la.gb.attribute(GATKVCFConstants.ULTIMA_INDEL_CLASSIFY, la.indel = (vc.getReference().length() < maxAlleleLength) ? C_INSERT : C_DELETE);

            /*
            lambda x: max([abs(len(y) - len(x['ref'])) for y in x['alleles']])
             */
            final int refLength = vc.getReference().length();
            la.gb.attribute(GATKVCFConstants.ULTIMA_INDEL_LENGTH, la.indelLength = vc.getAlleles().stream()
                    .filter(allele -> !allele.isReference())
                    .map(allele -> Math.abs(refLength - allele.length()))
                    .max(Integer::compare).get());
        }
    }

    // "hmer_indel_length" and "hmer_indel_nuc"
    private void isHmerIndel(final VariantContext vc, final LocalAttributes la) {

        if ( vc.isIndel() ) {
            List<Allele> alt = vc.getAlleles().stream()
                    .filter(allele -> !allele.isReference())
                    .collect(Collectors.toList());
            if ( alt.size() == 1 ) {
                Allele      altAllele = alt.get(0);
                if (C_INSERT.equals(la.indel)) {
                    /*
                    alt = [x for x in rec['alleles'] if x != rec['ref']][0][1:]
                    if len(set(alt)) != 1:
                        return (0, None)
                    elif fasta_idx[rec['chrom']][rec['pos']].seq.upper() != alt[0]:

                        return (0, None)
                    else:
                        return (utils.hmer_length(fasta_idx[rec['chrom']], rec['pos']), alt[0])
                     */
                    if (getReferenceNucleoid(la, vc.getStart()) == altAllele.getBases()[0]) {
                        la.gb.attribute(GATKVCFConstants.ULTIMA_HMER_INDEL_LENGTH, la.hmerIndelLength = hmerLength(la, vc.getStart()));
                        la.gb.attribute(GATKVCFConstants.ULTIMA_HMER_INDEL_NUC, Character.toString((char) altAllele.getBases()[0]));
                    }

                } else if (C_DELETE.equals(la.indel)) {
                    /*
                    del_seq = rec['ref'][1:]
                    if len(set(del_seq)) != 1:
                        return (0, None)
                    elif fasta_idx[rec['chrom']][rec['pos'] + len(rec['ref']) - 1].seq.upper() != del_seq[0]:
                        return (0, None)
                    else:
                        return (len(del_seq) + utils.hmer_length(fasta_idx[rec['chrom']],
                     */
                    if (getReferenceNucleoid(la,vc.getStart() + vc.getReference().length() - 1) == altAllele.getBases()[0]) {
                        la.gb.attribute(GATKVCFConstants.ULTIMA_HMER_INDEL_LENGTH, la.hmerIndelLength = altAllele.length() + hmerLength(la, vc.getStart() + vc.getReference().length() - 1));
                        la.gb.attribute(GATKVCFConstants.ULTIMA_HMER_INDEL_NUC, Character.toString((char) altAllele.getBases()[0]));
                    }
                }
            }
        }
    }

    private void getMotif(final VariantContext vc, final LocalAttributes la) {

        if ( la.indel != null && la.hmerIndelLength > 0 ) {
            annotateMotifAroundHherIndel(vc, la, la.hmerIndelLength);
        } else if ( la.indel != null && la.hmerIndelLength == 0 ) {
            annotateMotifAroundNonHherIndel(vc, la);
        } else {
            annotateMotifAroundSnp(vc, la);
        }
    }

    private void gcContent(final VariantContext vc, final LocalAttributes la) {

        /*
        chrom = faidx[rec['chrom']]
        beg = rec['pos'] - int(size / 2)
        end = beg + size
        seq = chrom[beg:end].seq.upper()
        seqGC = seq.replace('A', '').replace('T', '')
        return float(len(seqGC)) / len(seq)
         */
        int         begin = vc.getStart() - (GC_CONTENT_SIZE / 2);
        int         end = begin + GC_CONTENT_SIZE;
        String      seq = getReferenceMotif(la, begin, end);
        int         gcCount = 0;
        for ( byte b : seq.getBytes() ) {
            if ( b == 'G' || b == 'C' ) {
                gcCount++;
            }
        }
        la.gb.attribute(GATKVCFConstants.ULTIMA_GC_CONTENT, (float)gcCount / seq.length());
    }

    private void cycleSkip(final VariantContext vc, final LocalAttributes la) {

        /*
        #gt_field is None
        na_pos = df['indel'] | (df['alleles'].apply(len) > 2)
         */
        boolean naIndicator = la.indel != null || vc.getAlleles().size() > 2;

        /*
        df['cycleskip_status'] = "NA"
        snp_pos = ~na_pos
        snps = df.loc[snp_pos].copy()
        left_last = np.array(snps['left_motif']).astype(np.string_)
        right_first = np.array(snps['right_motif']).astype(np.string_)
         */
        String      css = C_CSS_NA;
        boolean     snpIndicator = !naIndicator;

        /*
        ref = np.array(snps['ref']).astype(np.string_)
        alt = np.array(snps['alleles'].apply(
            lambda x: x[1] if len(x) > 1 else None)).astype(np.string_)
         */
        String      ref = vc.getReference().getBaseString();
        String      alt = (snpIndicator && (vc.getAlleles().size() > 1)) ? vc.getAlleles().get(1).getBaseString() : null;

        /*
        ref_seqs = np.char.add(np.char.add(left_last, ref), right_first)
        alt_seqs = np.char.add(np.char.add(left_last, alt), right_first)
         */
        String      refSeqs = null;
        String      altSeqs = null;
        if ( la.leftMotif != null && la.rightMotif != null && alt != null ) {
            refSeqs = la.leftMotif + ref + la.rightMotif;
            altSeqs = la.leftMotif + alt + la.rightMotif;
        }

        /*
        ref_encs = [utils.generateKeyFromSequence(
                str(np.char.decode(x)), flow_order) for x in ref_seqs]
        alt_encs = [utils.generateKeyFromSequence(
                str(np.char.decode(x)), flow_order) for x in alt_seqs]
         */
        int[]   refEncs = generateKeyFromSequence(refSeqs, la.flowOrder);
        int[]   altEncs = generateKeyFromSequence(altSeqs, la.flowOrder);

        /*
        cycleskip = np.array([x for x in range(len(ref_encs))
                          if len(ref_encs[x]) != len(alt_encs[x])])
        poss_cycleskip = [x for x in range(len(ref_encs)) if len(ref_encs[x]) == len(alt_encs[x])
                          and (np.any(ref_encs[x][ref_encs[x] - alt_encs[x] != 0] == 0) or
                           np.any(alt_encs[x][ref_encs[x] - alt_encs[x] != 0] == 0))]
        s = set(np.concatenate((cycleskip, poss_cycleskip)))
        non_cycleskip = [x for x in range(len(ref_encs)) if x not in s]
         */
        if ( altEncs != null ) {
            boolean cycleskip = refEncs.length != altEncs.length;
            boolean passCycleskip = !cycleskip && !Arrays.equals(refEncs, altEncs);
            boolean s = cycleskip || passCycleskip;
            boolean nonCycleSkip = !s;

            if ( cycleskip )
                css = C_CSS_CS;
            else if ( passCycleskip )
                css = C_CSS_PCS;
            else if ( nonCycleSkip )
                css = C_CSS_NS;
        }

        la.gb.attribute(GATKVCFConstants.ULTIMA_CYCLESKIP_STATUS, css);
    }

    private int[] generateKeyFromSequence(final String sequence, final String flowOrder) {

        if ( sequence == null )
            return null;
        int[]         key = new int[sequence.length() * 4];
        /*
        key = []
        pos = 0
        for base in flow:
            hcount = 0
            for i in range(pos, len(sequence)):
                if sequence[i] == base:
                    hcount+=1
                else:
                    break
            else:
                key.append(hcount)
                break # end of sequence
            key.append(hcount)
            pos += hcount
         */
        byte[]      seq = sequence.getBytes();
        byte[]      flow = flowOrder.getBytes();
        int         pos = 0;
        int         keySize = 0;
        for ( int flowPos = 0 ; ; flowPos = (flowPos + 1) % flow.length ) {
            byte    base = flow[flowPos];
            int     hcount = 0;
            for ( int i = pos ; i < seq.length ; i++ ) {
                if ( seq[i] == 'N')
                    return null;
                else if ( seq[i] == base )
                    hcount++;
                else
                    break;
            }
            if ( pos >= seq.length ) {
                key[keySize++] = hcount;
                break;
            }
            key[keySize++] = hcount;
            pos += hcount;
        }

        return Arrays.copyOfRange(key, 0, keySize);
    }

    private void annotateMotifAroundHherIndel(final VariantContext vc, final LocalAttributes la, final int hmerLength) {
        /*
        chrom = faidx[rec['chrom']]
        pos = rec['pos']
        hmer_length = rec['hmer_indel_length']
        return chrom[pos - size:pos].seq.upper(), chrom[pos + hmer_length:pos + hmer_length + size].seq.upper()
         */
        la.gb.attribute(GATKVCFConstants.ULTIMA_LEFT_MOTIF, la.leftMotif = getReferenceMotif(la,vc.getStart() - MOTIF_SIZE, vc.getStart()));
        la.gb.attribute(GATKVCFConstants.ULTIMA_RIGHT_MOTIF, la.rightMotif = getReferenceMotif(la, vc.getStart() + hmerLength, vc.getStart() + hmerLength + MOTIF_SIZE));
    }

    private void annotateMotifAroundNonHherIndel(final VariantContext vc, final LocalAttributes la) {
        /*
        return chrom[pos - size:pos].seq.upper(),\
            chrom[pos + len(rec['ref']) - 1:pos +
                  len(rec['ref']) - 1 + size].seq.upper()
         */
        int         refLength = vc.getReference().length();
        la.gb.attribute(GATKVCFConstants.ULTIMA_LEFT_MOTIF, la.leftMotif = getReferenceMotif(la, vc.getStart() - MOTIF_SIZE, vc.getStart()));
        la.gb.attribute(GATKVCFConstants.ULTIMA_RIGHT_MOTIF, la.rightMotif = getReferenceMotif(la, vc.getStart() + refLength - 1, vc.getStart() + refLength - 1 + MOTIF_SIZE));
    }

    private void annotateMotifAroundSnp(final VariantContext vc, final LocalAttributes la) {
        /*
        return chrom[pos - size - 1:pos - 1].seq.upper(), chrom[pos:pos + size].seq.upper()
         */
        la.gb.attribute(GATKVCFConstants.ULTIMA_LEFT_MOTIF, la.leftMotif = getReferenceMotif(la, vc.getStart() - MOTIF_SIZE - 1, vc.getStart() - 1));
        la.gb.attribute(GATKVCFConstants.ULTIMA_RIGHT_MOTIF, la.rightMotif = getReferenceMotif(la, vc.getStart(), vc.getStart() + MOTIF_SIZE));
    }

    // get a single nucleoid from reference
    private byte getReferenceNucleoid(final LocalAttributes la, final int start) {
        int         index = start - la.ref.getWindow().getStart();
        byte[]      bases = la.ref.getBases();
        Utils.validIndex(index, bases.length);
        return bases[index];
    }

    // get motif from reference
    private String getReferenceMotif(final LocalAttributes la, final int start, final int end) {
        int         startIndex = start - la.ref.getWindow().getStart();
        int         endIndex = end - la.ref.getWindow().getStart();
        byte[]      bases = la.ref.getBases();
        Utils.validIndex(startIndex, bases.length);
        Utils.validIndex(endIndex, bases.length);
        return new String(Arrays.copyOfRange(bases, startIndex, endIndex));
    }

    private int hmerLength(final LocalAttributes la, final int start) {
        byte        base = getReferenceNucleoid(la, start);
        int         length = 1;
        while ( getReferenceNucleoid(la, start + length) == base )
            length++;
        return length;

    }

    @VisibleForTesting
    static Map<String, Object> annotateForTesting(final ReferenceContext ref, final VariantContext vc) {

        GenotypeBuilder             gb = new GenotypeBuilder();
        UltimaConcordanceAnnotator  annotator = new UltimaConcordanceAnnotator();

        annotator.annotate(ref, vc, null, gb, null);

        return gb.make().getExtendedAttributes();
    }

}
