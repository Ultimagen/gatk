package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.text.similarity.LevenshteinDistance;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.function.Consumer;

public abstract class BaseFeatureMapper {

    final int         surroundBefore;
    final int         surroundAfter;
    final int         minCigarElementLength;
    final LevenshteinDistance levDistance = new LevenshteinDistance();
    final Integer     smqSize;
    final Integer     smqSizeMean;

    final boolean     ignoreSurround;
    final int         spanBefore;
    final int         spanAfter;

    final FlowFeatureMapperArgumentCollection fmArgs;
    final SAMFileHeader hdr;

    BaseFeatureMapper(FlowFeatureMapperArgumentCollection fmArgs, SAMFileHeader hdr) {
        surroundBefore = fmArgs.snvIdenticalBases;
        surroundAfter = (fmArgs.snvIdenticalBasesAfter != 0) ?  fmArgs.snvIdenticalBasesAfter : surroundBefore;
        smqSize = fmArgs.surroundingMediaQualitySize;
        smqSizeMean = fmArgs.surroundingMeanQualitySize;
        this.fmArgs = fmArgs;
        this.hdr = hdr;

        // ignore surround
        ignoreSurround = fmArgs.reportAllAlts || fmArgs.tagBasesWithAdjacentRefDiff;
        spanBefore = ignoreSurround ? 0 : surroundBefore;
        spanAfter = ignoreSurround ? 0 : surroundAfter;
        minCigarElementLength = spanBefore + 1 + spanAfter;

        // adjust minimal read length
        FlowBasedRead.setMinimalReadLength(1 + 1 + spanAfter);
    }

    public FeatureMapper.FilterStatus noFeatureButFilterAt(GATKRead read, ReferenceContext referenceContext, int start) {

        // access bases
        final byte[]      bases = read.getBasesNoCopy();
        final byte[]      ref = referenceContext.getBases();

        // walk the cigar
        int         readOfs = 0;
        int         refOfs = 0;
        for ( final CigarElement cigarElement : read.getCigarElements() ) {

            final int     length = cigarElement.getLength();

            // worth looking into?
            boolean     includes = (start >= referenceContext.getStart() + refOfs) &&
                    (start < referenceContext.getStart() + refOfs + length);
            if ( includes && length >= minCigarElementLength &&
                    cigarElement.getOperator().consumesReadBases() &&
                    cigarElement.getOperator().consumesReferenceBases() ) {

                // break out if not enough clearing
                if ( (start < referenceContext.getStart() + refOfs + spanBefore) ||
                        (start >= referenceContext.getStart() + refOfs + length - spanAfter) )
                    return FeatureMapper.FilterStatus.Filtered;

                int         delta = start - (referenceContext.getStart() + refOfs);
                readOfs += delta;
                refOfs += delta;

                final boolean noFeature = bases[readOfs] == ref[refOfs];

                // check that this is really a SNV (must be surrounded by identical ref)
                boolean     surrounded = true;
                for ( int i = 0 ; i < surroundBefore && surrounded ; i++ ) {
                    final int bIndex = readOfs-1-i;
                    final int rIndex = refOfs-1-i;
                    if ( bIndex < 0 || bIndex >= bases.length || rIndex < 0 || rIndex >= ref.length ) {
                        surrounded = false;
                        continue;
                    }
                    if ( bases[bIndex] != ref[rIndex] ) {
                        surrounded = false;
                    }
                }
                for (int i = 0; i < surroundAfter && surrounded ; i++ ) {
                    final int bIndex = readOfs+1+i;
                    final int rIndex = refOfs+1+i;
                    if ( bIndex < 0 || bIndex >= bases.length || rIndex < 0 || rIndex >= ref.length ) {
                        surrounded = false;
                        continue;
                    }
                    if ( bases[bIndex] != ref[rIndex] ) {
                        surrounded = false;
                    }
                }
                if ( !surrounded ) {
                    continue;
                }

                return noFeature ? FeatureMapper.FilterStatus.NoFeatureAndFiltered : FeatureMapper.FilterStatus.Filtered;

            } else {

                // manual advance
                if (cigarElement.getOperator().consumesReadBases()) {
                    readOfs += length;
                }
                if (cigarElement.getOperator().consumesReferenceBases()) {
                    refOfs += length;
                }
            }
        };

        // if here, false
        return FeatureMapper.FilterStatus.None;
    }

    protected void reportFeatures(GATKRead read, ReferenceContext referenceContext, List<FlowFeatureMapper.MappedFeature> features, Consumer<? super FlowFeatureMapper.MappedFeature> action) {

        int refEditDistance = calcRefEditDistance(read, referenceContext);
        int nonIdentMBases = calcNonIndentBases(read, referenceContext);

        for ( FlowFeatureMapper.MappedFeature feature : features ) {
            feature.featuresOnRead = features.size();
            feature.nonIdentMBasesOnRead = nonIdentMBases;
            feature.refEditDistance = refEditDistance;
            action.accept(feature);
        }
    }

    protected void addSmq(FlowFeatureMapper.MappedFeature feature) {

        final GATKRead read = feature.read;

        // reverse quals?
        if ( smqSize != null || smqSizeMean != null  ) {

            // prepare qualities
            final byte[] quals;
            if (!read.isReverseStrand()) {
                quals = read.getBaseQualitiesNoCopy();
            } else {
                quals = read.getBaseQualities();
                ArrayUtils.reverse(quals);
            }

            // surrounding median quality?
            if ( smqSize != null ) {
                feature.smqLeft = calcSmq(quals, feature.index - 1 - smqSize, feature.index - 1, true);
                feature.smqRight = calcSmq(quals, feature.index + 1, feature.index + 1 + smqSize, true);
                if ( read.isReverseStrand() ) {

                    // left and right are reversed
                    int tmp = feature.smqLeft;
                    feature.smqLeft = feature.smqRight;
                    feature.smqRight = tmp;
                }
            }
            if ( smqSizeMean != null ) {
                feature.smqLeftMean = calcSmq(quals, feature.index - 1 - smqSizeMean, feature.index - 1, false);
                feature.smqRightMean = calcSmq(quals, feature.index + 1, feature.index + 1 + smqSizeMean, false);
                if ( read.isReverseStrand() ) {

                    // left and right are reversed
                    int tmp = feature.smqLeftMean;
                    feature.smqLeftMean = feature.smqRightMean;
                    feature.smqRightMean = tmp;
                }
            }
        }
    }

    protected int calcSmq(final byte[] quals, int from, int to, boolean median) {

        // limit from/to
        from = Math.max(0, Math.min(quals.length, from));
        to = Math.max(0, Math.min(quals.length, to - 1));
        if ( from > to ) {
            throw new GATKException("invalid qualities range: from > to");
        }

        // calc median
        byte[] range = Arrays.copyOfRange(quals, from, to + 1);
        if ( range.length == 0 ) {
            throw new GATKException("invalid qualities range: can't be empty");
        }

        if ( median ) {
            Arrays.sort(range);
            int midIndex = range.length / 2;
            if ((range.length % 2) == 1) {
                // odd
                return range[midIndex];
            } else {
                // even
                return (range[midIndex - 1] + range[midIndex]) / 2;
            }
        } else {
            int sum = 0;
            for ( int i = 0 ; i < range.length ; i++ ) {
                sum += range[i];
            }
            return sum / range.length;
        }
    }

    protected int calcRefEditDistance(GATKRead read, ReferenceContext referenceContext) {
        final byte[] bases = read.getBasesNoCopy();
        final byte[] ref = referenceContext.getBases();
        final int startSoftClip = read.getStart() - read.getSoftStart();
        final int endSoftClip = read.getSoftEnd() - read.getEnd();
        String basesString;
        if ( startSoftClip == 0 && endSoftClip == 0 ) {
            basesString = new String(bases);
        } else {
            basesString = new String(Arrays.copyOfRange(bases, startSoftClip, bases.length - endSoftClip));
        }
        return  levDistance.apply(basesString, new String(ref));
    }

    protected int calcNonIndentBases(GATKRead read, ReferenceContext referenceContext) {
        final byte[] bases = read.getBasesNoCopy();
        final byte[] ref = referenceContext.getBases();
        int nonIdentMBases = 0;
        int readOfs = 0;
        int refOfs = 0;
        for ( final CigarElement cigarElement : read.getCigarElements() ) {
            final int     length = cigarElement.getLength();
            if ( cigarElement.getOperator().consumesReadBases() && cigarElement.getOperator().consumesReferenceBases() ) {
                for ( int ofs = 0 ; ofs < length ; ofs++ ) {
                    if ( ref[refOfs+ofs] != 'N' && bases[readOfs+ofs] != ref[refOfs+ofs] ) {
                        nonIdentMBases++;
                    }
                }
            }
            if (cigarElement.getOperator().consumesReadBases()) {
                readOfs += length;
            }
            if (cigarElement.getOperator().consumesReferenceBases()) {
                refOfs += length;
            }
        }

        return nonIdentMBases;
    }

    protected boolean isSurrounded(GATKRead read, ReferenceContext referenceContext, final int readOfs, final int refOfs, int featureReadLength, int featureRefLength) {

        final byte[] bases = read.getBasesNoCopy();
        final byte[] ref = referenceContext.getBases();

        // check that this is really an isolated feature (must be surrounded by identical ref)
        boolean     surrounded = true;
        for ( int i = 0 ; i < surroundBefore && surrounded ; i++ ) {
            final int bIndex = readOfs-1-i;
            final int rIndex = refOfs-1-i;
            if ( bIndex < 0 || bIndex >= bases.length || rIndex < 0 || rIndex >= ref.length ) {
                surrounded = false;
                continue;
            }
            if ( bases[bIndex] != ref[rIndex] ) {
                surrounded = false;
            }
        }
        for (int i = 0; i < surroundAfter && surrounded ; i++ ) {
            final int bIndex = readOfs + featureReadLength + i;
            final int rIndex = refOfs + featureRefLength + i;
            if ( bIndex < 0 || bIndex >= bases.length || rIndex < 0 || rIndex >= ref.length ) {
                surrounded = false;
                continue;
            }
            if ( bases[bIndex] != ref[rIndex] ) {
                surrounded = false;
            }
        }

        return surrounded;
    }

    protected FlowFeatureMapper.MappedFeature enrichFeature(FlowFeatureMapper.MappedFeature feature, boolean surrounded) {
        final GATKRead read = feature.read;
        int hardLength = read.getUnclippedEnd() - read.getUnclippedStart() + 1;
        int readOfs = feature.readBasesOffset;
        feature.index = !read.isReverseStrand() ? readOfs : (hardLength - readOfs);
        if ( (fmArgs.reportAllAlts || fmArgs.tagBasesWithAdjacentRefDiff) && !surrounded )
            feature.adjacentRefDiff = true;
        addSmq(feature);
        return feature;
    }

    protected boolean ignoreBecauseNotSurrounded(final boolean surrounded) {
        return (!fmArgs.reportAllAlts && !fmArgs.tagBasesWithAdjacentRefDiff) && !surrounded;
    }

    public void forEachOnRead(GATKRead read, ReferenceContext referenceContext, Consumer<? super FlowFeatureMapper.MappedFeature> action) {

        // prepare list
        List<FlowFeatureMapper.MappedFeature>     features = new LinkedList<>();

        // walk the cigar (again) looking for features
        final byte[] bases = read.getBasesNoCopy();
        final byte[] ref = referenceContext.getBases();
        int readOfs = 0;
        int refOfs = 0;
        for ( final CigarElement cigarElement : read.getCigarElements() ) {

            final int     length = cigarElement.getLength();

            // worth looking into?
            if ( length >= minCigarElementLength &&
                    cigarElement.getOperator().consumesReadBases() &&
                    cigarElement.getOperator().consumesReferenceBases() ) {
                readOfs += spanBefore;
                refOfs += spanBefore;
                for ( int ofs = spanBefore ; ofs < length - spanAfter ; ofs++, readOfs++, refOfs++ ) {

                    if ( readOfs < bases.length && refOfs < ref.length ) {
                        FlowFeatureMapper.MappedFeature feature = detectFeature(read, referenceContext, bases, ref, readOfs, refOfs);
                        if (feature != null) {
                            features.add(feature);
                        }
                    }
                }
                readOfs += spanAfter;
                refOfs += spanAfter;

            } else {

                if ( isIndelMapper() ) {
                    FlowFeatureMapper.MappedFeature feature = null;
                    if ( cigarElement.getOperator() == CigarOperator.D ) {
                        feature = detectFeature(read, referenceContext, bases, ref, readOfs, refOfs, 0, length);
                    } else if ( cigarElement.getOperator() == CigarOperator.I ) {
                        feature = detectFeature(read, referenceContext, bases, ref, readOfs, refOfs, length, 0);
                    }
                    if ( feature != null ) {
                        features.add(feature);
                    }
                }

                // manual advance
                if (cigarElement.getOperator().consumesReadBases()) {
                    readOfs += length;
                }
                if (cigarElement.getOperator().consumesReferenceBases()) {
                    refOfs += length;
                }
            }
        };

        // report features
        reportFeatures(read, referenceContext, features, action);
    }

    protected FlowFeatureMapper.MappedFeature detectFeature(GATKRead read, ReferenceContext referenceContext, final byte bases[], final byte ref[], int readOfs, int refOfs) {
        return null;
    }

    protected FlowFeatureMapper.MappedFeature detectFeature(GATKRead read, ReferenceContext referenceContext, final byte bases[], final byte ref[], int readOfs, int refOfs, int readLength, int refLength) {
        return null;
    }
    protected boolean isIndelMapper() {
        return false;
    }
}