package uk.ac.ox.well.indiana.utils.io.cortex.graph;

public class CortexColor {
    private int meanReadLength;
    private long totalSequence;
    private String sampleName;
    private double errorRate;
    private boolean tipClippingApplied;
    private boolean lowCovgSupernodesRemoved;
    private boolean lowCovgKmersRemoved;
    private boolean cleanedAgainstGraph;
    private int lowCovSupernodesThreshold;
    private int lowCovKmerThreshold;
    private String cleanedAgainstGraphName;

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
