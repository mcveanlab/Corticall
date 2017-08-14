package uk.ac.ox.well.cortexjdk.utils.io.cortex.links;

import com.google.common.base.Joiner;

import java.util.*;

public class CortexJunctionsRecord {
    private boolean isFw;
    private int numKmers;
    private int numJunctions;
    private int[] coverages;
    private String junctions;
    private int[][] sources;

    public CortexJunctionsRecord(boolean isForwardOrientation, int numKmers, int numJunctions, int[] coverages, String junctions) {
        this.isFw = isForwardOrientation;
        this.numKmers = numKmers;
        this.numJunctions = numJunctions;
        this.coverages = coverages;
        this.junctions = junctions;
    }

    public CortexJunctionsRecord(boolean isForwardOrientation, int numKmers, int numJunctions, int[] coverages, String junctions, int[][] sources) {
        this.isFw = isForwardOrientation;
        this.numKmers = numKmers;
        this.numJunctions = numJunctions;
        this.coverages = coverages;
        this.junctions = junctions;
        this.sources = sources;
    }

    public String toString() {
        StringBuilder buffer = new StringBuilder();

        buffer.append(isFw ? "F" : "R").append(" ");
        buffer.append(numJunctions).append(" ");

        List<Integer> covs = new ArrayList<>();
        for (int coverage : coverages) {
            covs.add(coverage);
        }

        buffer.append(Joiner.on(",").join(covs)).append(" ");
        buffer.append(junctions);

        return buffer.toString();
    }

    public boolean isForward() { return isFw; }
    public boolean isFlipped() { return !isFw; }
    public int getNumKmers() { return numKmers; }
    public int getNumJunctions() { return numJunctions; }
    public int[] getCoverages() { return coverages; }
    public int getCoverage(int i) { return coverages[i]; }
    public String getJunctions() { return junctions; }
    public String getSources() {
        List<String> pieces = new ArrayList<>();

        for (int c = 0; c < sources.length; c++) {
            Set<Integer> sids = new TreeSet<>();

            for (int i = 0; i < sources[c].length; i++) {
                sids.add(sources[c][i]);
            }

            pieces.add(Joiner.on(",").join(sids));
        }

        return Joiner.on(";").join(pieces);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        CortexJunctionsRecord that = (CortexJunctionsRecord) o;

        if (isFw != that.isFw) return false;
        if (numJunctions != that.numJunctions) return false;
        if (numKmers != that.numKmers) return false;
        if (!Arrays.equals(coverages, that.coverages)) return false;
        if (!junctions.equals(that.junctions)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = (isFw ? 1 : 0);
        result = 31 * result + numKmers;
        result = 31 * result + numJunctions;
        result = 31 * result + Arrays.hashCode(coverages);
        result = 31 * result + junctions.hashCode();
        return result;
    }
}
