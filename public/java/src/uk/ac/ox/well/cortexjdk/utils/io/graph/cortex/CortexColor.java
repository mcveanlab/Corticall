package uk.ac.ox.well.cortexjdk.utils.io.graph.cortex;

public class CortexColor {
    private String sampleName;
    private int meanReadLength = 0;
    private long totalSequence = 0;
    private double errorRate = 0.0;
    private boolean tipClippingApplied = false;
    private boolean lowCovgSupernodesRemoved = false;
    private boolean lowCovgKmersRemoved = false;
    private boolean cleanedAgainstGraph = false;
    private int lowCovSupernodesThreshold = 0;
    private int lowCovKmerThreshold = 0;
    private String cleanedAgainstGraphName = "unknown";

    public int getMeanReadLength() {
        return meanReadLength;
    }

    public void setMeanReadLength(int meanReadLength) {
        this.meanReadLength = meanReadLength;
    }

    public long getTotalSequence() {
        return totalSequence;
    }

    public void setTotalSequence(long totalSequence) {
        this.totalSequence = totalSequence;
    }

    public String getSampleName() {
        return sampleName;
    }

    public void setSampleName(String sampleName) {
        this.sampleName = sampleName;
    }

    public double getErrorRate() {
        return errorRate;
    }

    public void setErrorRate(double errorRate) {
        this.errorRate = errorRate;
    }

    public boolean isTipClippingApplied() {
        return tipClippingApplied;
    }

    public void setTipClippingApplied(boolean tipClippingApplied) { this.tipClippingApplied = tipClippingApplied; }

    public boolean isLowCovgSupernodesRemoved() {
        return lowCovgSupernodesRemoved;
    }

    public void setLowCovgSupernodesRemoved(boolean lowCovgSupernodesRemoved) {
        this.lowCovgSupernodesRemoved = lowCovgSupernodesRemoved;
    }

    public boolean isLowCovgKmersRemoved() {
        return lowCovgKmersRemoved;
    }

    public void setLowCovgKmersRemoved(boolean lowCovgKmersRemoved) {
        this.lowCovgKmersRemoved = lowCovgKmersRemoved;
    }

    public boolean isCleanedAgainstGraph() {
        return cleanedAgainstGraph;
    }

    public void setCleanedAgainstGraph(boolean cleanedAgainstGraph) {
        this.cleanedAgainstGraph = cleanedAgainstGraph;
    }

    public int getLowCovSupernodesThreshold() {
        return lowCovSupernodesThreshold;
    }

    public void setLowCovSupernodesThreshold(int lowCovSupernodesThreshold) {
        this.lowCovSupernodesThreshold = lowCovSupernodesThreshold;
    }

    public int getLowCovKmerThreshold() {
        return lowCovKmerThreshold;
    }

    public void setLowCovKmerThreshold(int lowCovKmerThreshold) {
        this.lowCovKmerThreshold = lowCovKmerThreshold;
    }

    public String getCleanedAgainstGraphName() {
        return cleanedAgainstGraphName;
    }

    public void setCleanedAgainstGraphName(String cleanedAgainstGraphName) {
        this.cleanedAgainstGraphName = cleanedAgainstGraphName;
    }
}
