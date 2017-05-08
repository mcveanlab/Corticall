package uk.ac.ox.well.indiana.utils.traversal;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.stoppingconditions.TraversalStopper;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

public class TraversalEngineFactory {
    private TraversalEngineConfiguration configuration = new TraversalEngineConfiguration();

    public TraversalEngineFactory traversalColors(int... colors) { Arrays.stream(colors).forEach(c -> configuration.traversalColors.add(c)); return this; }
    public TraversalEngineFactory traversalColors(Collection<Integer> colors) { configuration.traversalColors.addAll(colors); return this; }
    public TraversalEngineFactory traversalSamples(String... samples) { Arrays.stream(samples).forEach(s -> configuration.traversalSamples.add(s)); return this; }
    public TraversalEngineFactory traversalSamples(Collection<String> samples) { configuration.traversalSamples.addAll(samples); return this; }

    public TraversalEngineFactory joiningColors(int... colors) { Arrays.stream(colors).forEach(c -> configuration.joiningColors.add(c)); return this; }
    public TraversalEngineFactory joiningColors(Collection<Integer> colors) { configuration.joiningColors.addAll(colors); return this; }
    public TraversalEngineFactory joiningSamples(String... samples) { Arrays.stream(samples).forEach(s -> configuration.joiningSamples.add(s)); return this; }
    public TraversalEngineFactory joiningSamples(Collection<String> samples) { configuration.joiningSamples.addAll(samples); return this; }
    public TraversalEngineFactory joiningGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> g) { return this; }

    public TraversalEngineFactory recruitmentColors(int... colors) { Arrays.stream(colors).forEach(c -> configuration.joiningColors.add(c)); return this; }
    public TraversalEngineFactory recruitmentColors(Collection<Integer> colors) { configuration.joiningColors.addAll(colors); return this; }
    public TraversalEngineFactory recruitmentSamples(String... samples) { Arrays.stream(samples).forEach(s -> configuration.recruitmentSamples.add(s)); return this; }
    public TraversalEngineFactory recruitmentSamples(Collection<String> samples) { configuration.recruitmentSamples.addAll(samples); return this; }

    public TraversalEngineFactory displayColors(int... colors) { Arrays.stream(colors).forEach(c -> configuration.joiningColors.add(c)); return this; }
    public TraversalEngineFactory displayColors(Collection<Integer> colors) { configuration.joiningColors.addAll(colors); return this; }
    public TraversalEngineFactory displaySamples(String... samples) { Arrays.stream(samples).forEach(s -> configuration.displaySamples.add(s)); return this; }
    public TraversalEngineFactory displaySamples(Collection<String> samples) { configuration.displaySamples.addAll(samples); return this; }

    public TraversalEngineFactory stopper(TraversalStopper<AnnotatedVertex, AnnotatedEdge> stoppingRule) { configuration.stoppingRule = stoppingRule; return this; }

    public TraversalEngineFactory links(CortexLinks... links) { Arrays.stream(links).forEach(l -> configuration.links.add(l)); return this; }
    public TraversalEngineFactory links(Collection<CortexLinks> links) { configuration.links.addAll(links); return this; }

    public TraversalEngineFactory graph(CortexGraph clean) { configuration.clean = clean; return this; }
    public TraversalEngineFactory dirty(CortexGraph dirty) { configuration.dirty = dirty; return this; }

    public TraversalEngine make() { return new TraversalEngine(configuration); }
}
