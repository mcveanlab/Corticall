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

    public GraphCombinationOperator gco = GraphCombinationOperator.OR;
    public TraversalDirection td = TraversalDirection.BOTH;
    public boolean connectUnusedNeighbors = false;

    public int traversalColor;
    public Set<Integer> joiningColors = new HashSet<>();
    public Set<Integer> recruitmentColors = new HashSet<>();
    public Set<Integer> displayColors = new HashSet<>();

    public DirectedGraph<CortexVertex, CortexEdge> previousGraph;
    public TraversalStopper<CortexVertex, CortexEdge> stoppingRule;

    public Map<Integer, CortexLinks> links = new HashMap<>();

    public CortexGraph graph;
    public CortexGraph dirty;
}
