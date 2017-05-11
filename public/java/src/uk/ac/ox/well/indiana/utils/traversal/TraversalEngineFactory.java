package uk.ac.ox.well.indiana.utils.traversal;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.stoppingconditions.TraversalStopper;

import java.util.Arrays;
import java.util.Collection;

public class TraversalEngineFactory {
    private TraversalEngineConfiguration configuration = new TraversalEngineConfiguration();

    public TraversalEngineFactory configuration(TraversalEngineConfiguration configuration) { this.configuration = configuration; return this; }

    public TraversalEngineFactory combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator op) { configuration.gco = op; return this; }
    public TraversalEngineFactory traversalDirection(TraversalEngineConfiguration.TraversalDirection td) { configuration.td = td; return this; }
    public TraversalEngineFactory connectUnusedNeighbors(boolean connectUnusedNeighbors) { configuration.connectUnusedNeighbors = connectUnusedNeighbors; return this; }

    public TraversalEngineFactory traversalColor(int color) { configuration.traversalColor = color; return this; }

    public TraversalEngineFactory joiningColors(int... colors) { Arrays.stream(colors).forEach(c -> configuration.joiningColors.add(c)); return this; }
    public TraversalEngineFactory joiningColors(Collection<Integer> colors) { configuration.joiningColors.addAll(colors); return this; }

    public TraversalEngineFactory previousGraph(DirectedGraph<CortexVertex, CortexEdge> previousGraph) { configuration.previousGraph = previousGraph; return this; }

    public TraversalEngineFactory recruitmentColors(int... colors) { Arrays.stream(colors).forEach(c -> configuration.recruitmentColors.add(c)); return this; }
    public TraversalEngineFactory recruitmentColors(Collection<Integer> colors) { configuration.recruitmentColors.addAll(colors); return this; }

    public TraversalEngineFactory displayColors(int... colors) { Arrays.stream(colors).forEach(c -> configuration.joiningColors.add(c)); return this; }
    public TraversalEngineFactory displayColors(Collection<Integer> colors) { configuration.joiningColors.addAll(colors); return this; }

    public TraversalEngineFactory stopper(TraversalStopper<CortexVertex, CortexEdge> stoppingRule) { configuration.stoppingRule = stoppingRule; return this; }

    //public TraversalEngineFactory links(CortexLinks... links) { Arrays.stream(links).forEach(l -> configuration.links.add(l)); return this; }
    //public TraversalEngineFactory links(Collection<CortexLinks> links) { configuration.links.addAll(links); return this; }

    public TraversalEngineFactory graph(CortexGraph clean) { configuration.graph = clean; return this; }
    //public TraversalEngineFactory dirty(CortexGraph dirty) { configuration.dirty = dirty; return this; }

    public TraversalEngine make() { return new TraversalEngine(configuration); }
}
