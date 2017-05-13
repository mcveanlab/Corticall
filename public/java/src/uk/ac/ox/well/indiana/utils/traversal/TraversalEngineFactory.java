package uk.ac.ox.well.indiana.utils.traversal;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.stoppingconditions.TraversalStopper;

import java.util.Arrays;
import java.util.Collection;

public class TraversalEngineFactory {
    private TraversalEngineConfiguration configuration = new TraversalEngineConfiguration();

    public TraversalEngineFactory configuration(TraversalEngineConfiguration configuration) { this.configuration = configuration; return this; }

    public TraversalEngineFactory combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator op) { configuration.setGraphCombinationOperator(op); return this; }
    public TraversalEngineFactory traversalDirection(TraversalEngineConfiguration.TraversalDirection td) { configuration.setTraversalDirection(td); return this; }
    public TraversalEngineFactory connectUnusedNeighbors(boolean connectUnusedNeighbors) { configuration.setConnectAllNeighbors(connectUnusedNeighbors); return this; }

    public TraversalEngineFactory traversalColor(int color) { configuration.setTraversalColor(color); return this; }

    public TraversalEngineFactory joiningColors(int... colors) { Arrays.stream(colors).forEach(c -> configuration.getJoiningColors().add(c)); return this; }
    public TraversalEngineFactory joiningColors(Collection<Integer> colors) { configuration.getJoiningColors().addAll(colors); return this; }

    public TraversalEngineFactory recruitmentColors(int... colors) { Arrays.stream(colors).forEach(c -> configuration.getRecruitmentColors().add(c)); return this; }
    public TraversalEngineFactory recruitmentColors(Collection<Integer> colors) { configuration.getRecruitmentColors().addAll(colors); return this; }

    public TraversalEngineFactory displayColors(int... colors) { Arrays.stream(colors).forEach(c -> configuration.getJoiningColors().add(c)); return this; }
    public TraversalEngineFactory displayColors(Collection<Integer> colors) { configuration.getJoiningColors().addAll(colors); return this; }

    public TraversalEngineFactory previousTraversal(DirectedGraph<CortexVertex, CortexEdge> previousTraversal) { configuration.setPreviousTraversal(previousTraversal); return this; }

    public TraversalEngineFactory stopper(TraversalStopper<CortexVertex, CortexEdge> stoppingRule) { configuration.setStoppingRule(stoppingRule); return this; }

    //public TraversalEngineFactory links(CortexLinks... links) { Arrays.stream(links).forEach(l -> configuration.links.add(l)); return this; }
    //public TraversalEngineFactory links(Collection<CortexLinks> links) { configuration.links.addAll(links); return this; }

    public TraversalEngineFactory graph(CortexGraph clean) { configuration.setGraph(clean); return this; }
    //public TraversalEngineFactory dirty(CortexGraph dirty) { configuration.dirty = dirty; return this; }

    public TraversalEngine make() { return new TraversalEngine(configuration); }
}
