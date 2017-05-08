package uk.ac.ox.well.indiana.utils.traversal;

import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.stoppingconditions.TraversalStopper;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by kiran on 05/05/2017.
 */
public class TraversalEngineConfiguration {
    public Set<Integer> traversalColors = new HashSet<>();
    public Set<Integer> joiningColors = new HashSet<>();
    public Set<Integer> recruitmentColors = new HashSet<>();
    public Set<Integer> displayColors = new HashSet<>();

    public Set<String> traversalSamples = new HashSet<>();
    public Set<String> joiningSamples = new HashSet<>();
    public Set<String> recruitmentSamples = new HashSet<>();
    public Set<String> displaySamples = new HashSet<>();

    public Set<CortexLinks> links = new HashSet<>();

    public CortexGraph clean;
    public CortexGraph dirty;

    public TraversalStopper<AnnotatedVertex, AnnotatedEdge> stoppingRule;

    public boolean useDirtyGraph;
}
