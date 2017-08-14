package uk.ac.ox.well.cortexjdk.utils.traversal;

import org.jgrapht.DirectedGraph;
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
    private Set<Integer> displayColors = new TreeSet<>();

    private DirectedGraph<CortexVertex, CortexEdge> previousTraversal;
    private Class<? extends TraversalStopper<CortexVertex, CortexEdge>> stoppingRule;

    //private Map<Integer, CortexLinks> links = new HashMap<>();
    //private Map<CortexLinksMap, String> links = new HashMap<>();
    private Map<CortexLinks, String> links = new HashMap<>();

    private CortexGraph graph;
    private CortexGraph rois;

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

    public Set<Integer> getDisplayColors() { return displayColors; }
    public void setDisplayColors(Collection<Integer> displayColors) { this.displayColors = new TreeSet<>(displayColors); }
    public void setDisplayColors() { this.displayColors.clear(); }

    public DirectedGraph<CortexVertex, CortexEdge> getPreviousTraversal() { return previousTraversal; }
    public void setPreviousTraversal(DirectedGraph<CortexVertex, CortexEdge> previousTraversal) { this.previousTraversal = previousTraversal; }

    public Class<? extends TraversalStopper<CortexVertex, CortexEdge>> getStoppingRule() { return stoppingRule; }
    public void setStoppingRule(Class<? extends TraversalStopper<CortexVertex, CortexEdge>> stoppingRule) { this.stoppingRule = stoppingRule; }

    //public Map<Integer, CortexLinks> getLinks() { return links; }
    //public void setLinks(Map<Integer, CortexLinks> links) { this.links = links; }

    //public Map<CortexLinksMap, String> getLinks() { return links; }
    //public void setLinks(Map<CortexLinksMap, String> links) { this.links = links; }

    public Map<CortexLinks, String> getLinks() { return links; }
    public void setLinks(Map<CortexLinks, String> links) { this.links = links; }

    public CortexGraph getGraph() { return graph; }
    public void setGraph(CortexGraph graph) { this.graph = graph; }
    //public void setGraph(CortexGraph graph) { this.graph = new CortexGraph(graph.getCortexFile()); }

    public CortexGraph getRois() { return rois; }
    public void setRois(CortexGraph rois) { this.rois = rois; }
    //public void setRois(CortexGraph rois) { this.rois = new CortexGraph(rois.getCortexFile()); }
}
