package uk.ac.ox.well.cortexjdk.utils.caller;

import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinks;

import java.util.*;

/**
 * Created by kiran on 30/08/2017.
 */
public class BubbleCallerConfiguration {
    private int alternateColor = -1;
    private Set<Integer> referenceColors = new TreeSet<>();

    private DeBruijnGraph graph;
    private DeBruijnGraph rois;
    private Set<CortexLinks> links = new HashSet<>();
    private Map<String, KmerLookup> references;

    public int getAlternateColor() { return alternateColor; }
    public void setAlternateColor(int alternateColor) { this.alternateColor = alternateColor; }
    public void setAlternateColor() { this.alternateColor = -1; }

    public Set<Integer> getReferenceColors() { return referenceColors; }
    public void setReferenceColors(Collection<Integer> referenceColors) { this.referenceColors = new TreeSet<>(referenceColors); }
    public void setReferenceColors() { this.referenceColors.clear(); }

    public Set<CortexLinks> getLinks() { return links; }
    public void setLinks(Set<CortexLinks> links) { this.links = links; }

    public DeBruijnGraph getGraph() { return graph; }
    public void setGraph(DeBruijnGraph graph) { this.graph = graph; }

    public DeBruijnGraph getRois() { return rois; }
    public void setRois(DeBruijnGraph rois) { this.rois = rois; }

    public Map<String, KmerLookup> getReferences() { return this.references; }
    public void setReferences(Map<String, KmerLookup> references) { this.references = references; }
}
