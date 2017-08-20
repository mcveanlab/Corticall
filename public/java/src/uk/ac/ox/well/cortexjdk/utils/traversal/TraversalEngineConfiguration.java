package uk.ac.ox.well.cortexjdk.utils.traversal;

import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.stoppingconditions.TraversalStopper;

import java.util.*;

/**
 * Created by kiran on 05/05/2017.
 */
public class TraversalEngineConfiguration {
    public enum GraphCombinationOperator { AND, OR }
    public enum TraversalDirection { BOTH, FORWARD, REVERSE }

    private GraphCombinationOperator gco = GraphCombinationOperator.AND;
    private TraversalDirection td = TraversalDirection.BOTH;
    private boolean connectAllNeighbors = false;

    private int traversalColor = -1;
    private Set<Integer> joiningColors = new TreeSet<>();
    private Set<Integer> recruitmentColors = new TreeSet<>();
    private Set<Integer> secondaryColors = new TreeSet<>();

    private DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousTraversal;
    private Class<? extends TraversalStopper<CortexVertex, CortexEdge>> stoppingRule;

    private DeBruijnGraph graph;
    private DeBruijnGraph rois;
    private Set<CortexLinks> links = new HashSet<>();

    public GraphCombinationOperator getGraphCombinationOperator() { return gco; }
    public void setGraphCombinationOperator(GraphCombinationOperator gco) { this.gco = gco; }

    public TraversalDirection getTraversalDirection() { return td; }
    public void setTraversalDirection(TraversalDirection td) { this.td = td; }

    public boolean connectAllNeighbors() { return connectAllNeighbors; }
    public void setConnectAllNeighbors(boolean connectAllNeighbors) { this.connectAllNeighbors = connectAllNeighbors; }

    public int getTraversalColor() { return traversalColor; }
    public void setTraversalColor(int traversalColor) { this.traversalColor = traversalColor; }
    public void setTraversalColor() { this.traversalColor = -1; }

    public Set<Integer> getJoiningColors() { return joiningColors; }
    public void setJoiningColors(Collection<Integer> joiningColors) { this.joiningColors = new TreeSet<>(joiningColors); }
    public void setJoiningColors() { this.joiningColors.clear(); }

    public Set<Integer> getRecruitmentColors() { return recruitmentColors; }
    public void setRecruitmentColors(Collection<Integer> recruitmentColors) { this.recruitmentColors = new TreeSet<>(recruitmentColors); }
    public void setRecruitmentColors() { this.recruitmentColors.clear(); }

    public Set<Integer> getSecondaryColors() { return secondaryColors; }
    public void setSecondaryColors(Collection<Integer> secondaryColors) { this.secondaryColors = new TreeSet<>(secondaryColors); }
    public void setSecondaryColors() { this.secondaryColors.clear(); }

    public DirectedWeightedPseudograph<CortexVertex, CortexEdge> getPreviousTraversal() { return previousTraversal; }
    public void setPreviousTraversal(DirectedWeightedPseudograph<CortexVertex, CortexEdge> previousTraversal) { this.previousTraversal = previousTraversal; }

    public Class<? extends TraversalStopper<CortexVertex, CortexEdge>> getStoppingRule() { return stoppingRule; }
    public void setStoppingRule(Class<? extends TraversalStopper<CortexVertex, CortexEdge>> stoppingRule) { this.stoppingRule = stoppingRule; }

    public Set<CortexLinks> getLinks() { return links; }
    public void setLinks(Set<CortexLinks> links) { this.links = links; }

    public DeBruijnGraph getGraph() { return graph; }
    public void setGraph(DeBruijnGraph graph) { this.graph = graph; }

    public DeBruijnGraph getRois() { return rois; }
    public void setRois(DeBruijnGraph rois) { this.rois = rois; }
}
