package uk.ac.ox.well.indiana.utils.traversal;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.stoppingconditions.TraversalStopper;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created by kiran on 05/05/2017.
 */
public class TraversalEngineConfiguration {
    public enum GraphCombinationOperator { OR, AND }
    public enum TraversalDirection { BOTH, FORWARD, REVERSE }

    private GraphCombinationOperator gco = GraphCombinationOperator.OR;
    private TraversalDirection td = TraversalDirection.BOTH;
    private boolean connectAllNeighbors = false;

    private int traversalColor = -1;
    private Set<Integer> joiningColors = new HashSet<>();
    private Set<Integer> recruitmentColors = new HashSet<>();
    private Set<Integer> displayColors = new HashSet<>();

    private DirectedGraph<CortexVertex, CortexEdge> previousTraversal;
    private TraversalStopper<CortexVertex, CortexEdge> stoppingRule;

    private Map<Integer, CortexLinks> links = new HashMap<>();

    private CortexGraph graph;
    private CortexGraph dirty;

    public GraphCombinationOperator getGraphCombinationOperator() {
        return gco;
    }

    public void setGraphCombinationOperator(GraphCombinationOperator gco) {
        this.gco = gco;
    }

    public TraversalDirection getTraversalDirection() {
        return td;
    }

    public void setTraversalDirection(TraversalDirection td) {
        this.td = td;
    }

    public boolean connectAllNeighbors() {
        return connectAllNeighbors;
    }

    public void setConnectAllNeighbors(boolean connectAllNeighbors) {
        this.connectAllNeighbors = connectAllNeighbors;
    }

    public int getTraversalColor() {
        return traversalColor;
    }

    public void setTraversalColor(int traversalColor) {
        this.traversalColor = traversalColor;
    }

    public Set<Integer> getJoiningColors() {
        return joiningColors;
    }

    public void setJoiningColors(Set<Integer> joiningColors) {
        this.joiningColors = joiningColors;
    }

    public Set<Integer> getRecruitmentColors() {
        return recruitmentColors;
    }

    public void setRecruitmentColors(Set<Integer> recruitmentColors) {
        this.recruitmentColors = recruitmentColors;
    }

    public Set<Integer> getDisplayColors() {
        return displayColors;
    }

    public void setDisplayColors(Set<Integer> displayColors) {
        this.displayColors = displayColors;
    }

    public DirectedGraph<CortexVertex, CortexEdge> getPreviousTraversal() {
        return previousTraversal;
    }

    public void setPreviousTraversal(DirectedGraph<CortexVertex, CortexEdge> previousTraversal) {
        this.previousTraversal = previousTraversal;
    }

    public TraversalStopper<CortexVertex, CortexEdge> getStoppingRule() {
        return stoppingRule;
    }

    public void setStoppingRule(TraversalStopper<CortexVertex, CortexEdge> stoppingRule) {
        this.stoppingRule = stoppingRule;
    }

    public Map<Integer, CortexLinks> getLinks() {
        return links;
    }

    public void setLinks(Map<Integer, CortexLinks> links) {
        this.links = links;
    }

    public CortexGraph getGraph() {
        return graph;
    }

    public void setGraph(CortexGraph graph) {
        this.graph = graph;
    }

    public CortexGraph getDirty() {
        return dirty;
    }

    public void setDirty(CortexGraph dirty) {
        this.dirty = dirty;
    }
}
