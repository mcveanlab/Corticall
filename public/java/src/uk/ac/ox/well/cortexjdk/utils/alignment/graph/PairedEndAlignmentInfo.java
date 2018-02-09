package uk.ac.ox.well.cortexjdk.utils.alignment.graph;

import java.util.Map;
import java.util.Set;

public class PairedEndAlignmentInfo {
    private SingleEndAlignmentInfo sai1;
    private SingleEndAlignmentInfo sai2;
    private Map<Integer, Set<Integer>> dists;

    public PairedEndAlignmentInfo(SingleEndAlignmentInfo sai1, SingleEndAlignmentInfo sai2, Map<Integer, Set<Integer>> dists) {
        this.sai1 = sai1;
        this.sai2 = sai2;
        this.dists = dists;
    }

    public SingleEndAlignmentInfo getFirstEndAlignment() { return sai1; }

    public SingleEndAlignmentInfo getSecondEndAlignment() { return sai2; }

    public Map<Integer, Set<Integer>> getAlignmentDistanceMap() { return dists; }
}
