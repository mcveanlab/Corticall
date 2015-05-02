package uk.ac.ox.well.indiana.utils.io.cortex.links;

import com.google.common.base.Joiner;

import java.util.Arrays;

public class CortexJunctionsRecord {
    private boolean isFw;
    private int numKmers;
    private int numJunctions;
    private int[] coverages;
    private String junctions;
    private String seq;
    private int[] juncpos;

    public CortexJunctionsRecord(boolean isForwardOrientation, int numKmers, int numJunctions, int[] coverages, String junctions, String seq, int[] juncpos) {
        this.isFw = isForwardOrientation;
        this.numKmers = numKmers;
        this.numJunctions = numJunctions;
        this.coverages = coverages;
        this.junctions = junctions;
        this.seq = seq;
        this.juncpos = juncpos;
    }

    public String toString() {
        StringBuilder buffer = new StringBuilder();

        buffer.append(isFw ? "F" : "R").append(" ");
        buffer.append(numKmers).append(" ");
        buffer.append(numJunctions).append(" ");

        for (int c = 0; c < coverages.length; c++) {
            buffer.append(coverages[c]).append(" ");
        }

        buffer.append(junctions);

        if (seq != null) {
            buffer.append(" seq=").append(seq);
        }

        if (juncpos != null) {
            buffer.append(" juncpos=");

            for (int j = 0; j < juncpos.length; j++) {
                buffer.append(juncpos[j]);

                if (j < juncpos.length - 1) {
                    buffer.append(",");
                }
            }
        }

        return buffer.toString();
    }

    public boolean isForward() { return isFw; }
    public boolean isFlipped() { return !isFw; }
    public int getNumKmers() { return numKmers; }
    public int getNumJunctions() { return numJunctions; }
    public int[] getCoverages() { return coverages; }
    public int getCoverage(int i) { return coverages[i]; }
    public String getJunctions() { return junctions; }
    public String getSeq() { return seq; }
    public int[] getJunctionPositions() { return juncpos; }

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
