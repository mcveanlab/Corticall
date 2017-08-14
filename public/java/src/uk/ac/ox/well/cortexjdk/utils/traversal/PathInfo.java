package uk.ac.ox.well.cortexjdk.utils.traversal;

public class PathInfo {
    public String start, stop, child, parent;

    public PathInfo(String start, String stop, String child, String parent) {
        this.start = start;
        this.stop = stop;
        this.child = child;
        this.parent = parent;
    }
}
