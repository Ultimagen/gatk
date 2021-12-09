package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;

public class Feature implements Comparable<Feature> {

    GATKRead    read;
    FlowFeatureMapperArgumentCollection.MappingFeatureEnum  type;
    byte[]      readBases;
    byte[]      refBases;
    int         readBasesOffset; // offset of read bases array
    int         start;      // location (on rerence)
    int         offsetDelta;
    double      score;
    int         readCount;
    int         filteredCount;
    int         nonIdentMBasesOnRead;
    int         featuresOnRead;
    int         refEditDistance;
    int         index;

    public Feature(GATKRead read, FlowFeatureMapperArgumentCollection.MappingFeatureEnum  type, byte[] readBases,
                   byte[] refBases, int readBasesOffset, int start, int offsetDelta) {
        this.read = read;
        this.type = type;
        this.readBases = readBases;
        this.refBases = refBases;
        this.readBasesOffset = readBasesOffset;
        this.start = start;
        this.offsetDelta = offsetDelta;
    }

    static Feature makeSNV(GATKRead read, int offset, byte refBase, int start, int offsetDelta) {
        byte[]      readBases = {read.getBasesNoCopy()[offset]};
        byte[]      refBases = {refBase};
        return new Feature(
                read,
                FlowFeatureMapperArgumentCollection.MappingFeatureEnum.SNV,
                readBases,
                refBases,
                offset,
                start,
                offsetDelta);
    }

    @Override
    public String toString() {
        return "Feature{" +
                "read=" + read +
                ", type=" + type +
                ", readBases=" + Arrays.toString(readBases) +
                ", refBases=" + Arrays.toString(refBases) +
                ", readBasesOffset=" + readBasesOffset +
                ", start=" + start +
                '}';
    }

    @Override
    public int compareTo(Feature o) {

        int     delta = this.read.getContig().compareTo(o.read.getContig());
        if ( delta != 0 )
            return delta;

        delta = this.start - o.start;

        return delta;
    }
}
