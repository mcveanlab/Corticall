package uk.ac.ox.well.cortexjdk.utils.traversal;

import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.io.graph.ConnectivityAnnotations;
import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.TraversalStoppingRule;

import java.util.*;
import java.util.stream.Stream;

/**
 * Created by kiran on 05/05/2017.
 */
public class TraversalEngineConfiguration {
    public enum GraphCombinationOperator { AND, OR }
    public enum TraversalDirection { BOTH, FORWARD, REVERSE }

    private Set<Integer> traversalColors = new LinkedHashSet<>();
    private Set<Integer> joiningColors = new TreeSet<>();
    private Set<Integer> recruitmentColors = new TreeSet<>();
    private Set<Integer> secondaryColors = new TreeSet<>();

    private GraphCombinationOperator gco = GraphCombinationOperator.OR;
    private TraversalDirection td = TraversalDirection.BOTH;
    private boolean connectAllNeighbors = false;

    private int maxLength = 75000;

    private Class<? extends TraversalStoppingRule<CortexVertex, CortexEdge>> stoppingRule = ContigStopper.class;

    private DeBruijnGraph graph;
    private DeBruijnGraph rois;
    private Set<ConnectivityAnnotations> links = new HashSet<>();
    private Set<IndexedReference> kls = new HashSet<>();

    private boolean debug = false;

    public GraphCombinationOperator getGraphCombinationOperator() { return gco; }
    public void setGraphCombinationOperator(GraphCombinationOperator gco) { this.gco = gco; }

    public TraversalDirection getTraversalDirection() { return td; }
    public void setTraversalDirection(TraversalDirection td) { this.td = td; }

    public boolean connectAllNeighbors() { return connectAllNeighbors; }
    public void setConnectAllNeighbors(boolean connectAllNeighbors) { this.connectAllNeighbors = connectAllNeighbors; }

    public Set<Integer> getTraversalColors() { return traversalColors; }
    public void setTraversalColors(Collection<Integer> traversalColors) { this.traversalColors = new LinkedHashSet<>(traversalColors); }
    public void setTraversalColors(int... traversalColors) { Arrays.stream(traversalColors).forEach(c -> this.traversalColors.add(c)); }

    public Set<Integer> getJoiningColors() { return joiningColors; }
    public void setJoiningColors(Collection<Integer> joiningColors) { this.joiningColors = new TreeSet<>(joiningColors); }
    public void setJoiningColors(int... joiningColors) { Arrays.stream(joiningColors).forEach(c -> this.joiningColors.add(c)); }

    public Set<Integer> getRecruitmentColors() { return recruitmentColors; }
    public void setRecruitmentColors(Collection<Integer> recruitmentColors) { this.recruitmentColors = new TreeSet<>(recruitmentColors); }
    public void setRecruitmentColors(int... recruitmentColors) { Arrays.stream(recruitmentColors).forEach(c -> this.recruitmentColors.add(c)); }

    public Set<Integer> getSecondaryColors() { return secondaryColors; }
    public void setSecondaryColors(Collection<Integer> secondaryColors) { this.secondaryColors = new TreeSet<>(secondaryColors); }
    public void setSecondaryColors(int... secondaryColors) { Arrays.stream(secondaryColors).forEach(c -> this.secondaryColors.add(c)); }

    public Class<? extends TraversalStoppingRule<CortexVertex, CortexEdge>> getStoppingRule() { return stoppingRule; }
    public void setStoppingRule(Class<? extends TraversalStoppingRule<CortexVertex, CortexEdge>> stoppingRule) { this.stoppingRule = stoppingRule; }

    public Set<ConnectivityAnnotations> getLinks() { return links; }
    public void setLinks(Set<ConnectivityAnnotations> links) { this.links = links; }

    public DeBruijnGraph getGraph() { return graph; }
    public void setGraph(DeBruijnGraph graph) { this.graph = graph; }

    public DeBruijnGraph getRois() { return rois; }
    public void setRois(DeBruijnGraph rois) { this.rois = rois; }

    public Set<IndexedReference> getReferences() { return kls; }
    public void setReferences(Set<IndexedReference> kls) { this.kls = kls; }

    public void setMaxWalkLength(int maxLength) { this.maxLength = maxLength; }
    public int getMaxBranchLength() { return maxLength; }

    public void setDebugFlag() { this.debug = true; }
    public boolean getDebugFlag() { return this.debug; }
}
