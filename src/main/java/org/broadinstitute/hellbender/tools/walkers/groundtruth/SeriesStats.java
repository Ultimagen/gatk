package org.broadinstitute.hellbender.tools.walkers.groundtruth;

import org.apache.commons.collections.map.LazySortedMap;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicInteger;

public class SeriesStats {

    // local state
    private double last = Double.NaN;
    private int count = 0;
    private double sum = 0;
    private double min = Double.NaN;
    private double max = Double.NaN;
    private SortedMap<Double, AtomicInteger> bins = new TreeMap<>();
    private int intCount = 0;
    private Map<Double, SeriesStats> auxBins = new LinkedHashMap<>();

    public void csvWrite(final String path) throws IOException {
        PrintWriter pw = new PrintWriter(path);
        pw.println("value,count");
        boolean intKeys = (count == intCount);
        for (Map.Entry<Double, AtomicInteger> entry : bins.entrySet() ) {
            if ( intKeys ) {
                pw.println(String.format("%d,%d", entry.getKey().intValue(), entry.getValue().get()));
            } else {
                pw.println(String.format("%f,%d", entry.getKey(), entry.getValue().get()));
            }
        }
        pw.close();
    }

    public void add(int v) {
        add((double)v);
        intCount++;
    }

    public void add(double v) {

        // save in simple values
        last = v;
        sum += v;
        if ( count > 0 ) {
            min = Math.min(min, v);
            max = Math.max(max, v);
        } else {
            min = max = v;
        }
        count++;

        // save in bins
        if ( bins.containsKey(v) ) {
            bins.get(v).incrementAndGet();
        } else {
            bins.put(v, new AtomicInteger(1));
        }
    }

    public void aux(int v, int auxValue) {
        aux((double)v, auxValue);
    }

    public void aux(double v, int auxValue) {

        if ( auxBins.containsKey(v) ) {
            auxBins.get(v).add(auxValue);
        } else {
            SeriesStats ss = new SeriesStats();
            ss.add(auxValue);
            auxBins.put(v, ss);
        }
    }

    public double getLast() {
        return last;
    }

    public int getCount() {
        return count;
    }

    public double getMin() {
        return (count != 0) ? min : Double.NaN;
    }

    public double getMax() {
        return (count != 0) ? max : Double.NaN;
    }

    public int getUniq() {
        return bins.size();
    }

    public double getMean() {
        return (count != 0) ? (sum / count) : Double.NaN;
    }

    public double getMedian() {
        return getPercentile(50);
    }

    public double getPercentile(double precentile) {
        if ( count == 0 ) {
            return Double.NaN;
        } else if ( count == 1 ) {
            return last;
        } else {

            int percentileIndex = (int)(count * precentile / 100);
            int index = 0;
            for (Map.Entry<Double, AtomicInteger> entry : bins.entrySet() ) {
                int binSize = entry.getValue().get();
                if ( percentileIndex >= index && (percentileIndex < (index + binSize)) ) {
                    return entry.getKey();
                }
                index += binSize;
            }

            // if here, we need the highest entry
            return bins.lastKey();
        }
    }

    public double getStd() {

        if (count == 0) {
            return Double.NaN;
        }

        // calculate mean
        double mean = getMean();

        // calculate sum of sqr(deviations)
        double vsum = 0;
        for (Map.Entry<Double, AtomicInteger> entry : bins.entrySet()) {
            int binSize = entry.getValue().get();
            vsum += (Math.pow(entry.getKey() - mean, 2) * binSize);
        }

        // calculate variance and std
        double variance = vsum / count;
        return Math.sqrt(variance);
    }

    public Map<Double, AtomicInteger> getBins() {
        return this.bins;
    }

    public Map<Double, SeriesStats> getAuxBins() {
        return this.auxBins;
    }

    public String toDigest() {
        return this.bins.toString();
    }
}