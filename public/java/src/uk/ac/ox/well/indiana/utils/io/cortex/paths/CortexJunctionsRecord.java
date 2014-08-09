package uk.ac.ox.well.indiana.utils.io.cortex.paths;

import java.util.Arrays;

public class CortexJunctionsRecord {
    private boolean isFw;
    private int numKmers;
    private int numJunctions;
    private int[] coverages;
    private String junctions;

    public CortexJunctionsRecord(boolean isForwardOrientation, int numKmers, int numJunctions, int[] coverages, String junctions) {
        this.isFw = isForwardOrientation;
        this.numKmers = numKmers;
        this.numJunctions = numJunctions;
        this.coverages = coverages;
        this.junctions = junctions;
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

        return buffer.toString();
    }

    public boolean isForwardOrientation() { return isFw; }
    public int getNumKmers() { return numKmers; }
    public int getNumJunctions() { return numJunctions; }
    public int[] getCoverages() { return coverages; }
    public int getCoverage(int i) { return coverages[i]; }
    public String getJunctions() { return junctions; }

    public boolean equals(Object obj) {
        if (obj instanceof CortexJunctionsRecord) {
            return (isFw == ((CortexJunctionsRecord) obj).isForwardOrientation()) &&
                   (numKmers == ((CortexJunctionsRecord) obj).numKmers) &&
                   (numJunctions == ((CortexJunctionsRecord) obj).numJunctions) &&
                   (junctions.equals(((CortexJunctionsRecord) obj).getJunctions())) &&
                   (Arrays.equals(coverages, ((CortexJunctionsRecord) obj).getCoverages()));
        }

        return false;
    }
}
