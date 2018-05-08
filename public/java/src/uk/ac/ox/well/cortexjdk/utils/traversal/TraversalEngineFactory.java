package uk.ac.ox.well.cortexjdk.utils.traversal;

import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.TraversalStoppingRule;

import java.util.Arrays;
import java.util.Collection;

public class TraversalEngineFactory {
    private TraversalEngineConfiguration configuration = new TraversalEngineConfiguration();

    public TraversalEngineFactory configuration(TraversalEngineConfiguration configuration) { this.configuration = configuration; return this; }

    public TraversalEngineFactory combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator op) { configuration.setGraphCombinationOperator(op); return this; }
    public TraversalEngineFactory traversalDirection(TraversalEngineConfiguration.TraversalDirection td) { configuration.setTraversalDirection(td); return this; }

    public TraversalEngineFactory connectAllNeighbors(boolean connectAllNeighbors) { configuration.setConnectAllNeighbors(connectAllNeighbors); return this; }
    public TraversalEngineFactory maxBranchLength(int maxLength) { configuration.setMaxWalkLength(maxLength); return this; }

    public TraversalEngineFactory traversalColors() { configuration.getTraversalColors().clear(); return this; }
    public TraversalEngineFactory traversalColors(int... colors) { Arrays.stream(colors).forEach(c -> configuration.getTraversalColors().add(c)); return this; }
    public TraversalEngineFactory traversalColors(Collection<Integer> colors) { configuration.getTraversalColors().addAll(colors); return this; }

    public TraversalEngineFactory joiningColors() { configuration.getJoiningColors().clear(); return this; }
    public TraversalEngineFactory joiningColors(int... colors) { Arrays.stream(colors).forEach(c -> configuration.getJoiningColors().add(c)); return this; }
    public TraversalEngineFactory joiningColors(Collection<Integer> colors) { configuration.getJoiningColors().addAll(colors); return this; }

    public TraversalEngineFactory recruitmentColors() { configuration.getRecruitmentColors().clear(); return this; }
    public TraversalEngineFactory recruitmentColors(int... colors) { Arrays.stream(colors).forEach(c -> configuration.getRecruitmentColors().add(c)); return this; }
    public TraversalEngineFactory recruitmentColors(Collection<Integer> colors) { configuration.getRecruitmentColors().addAll(colors); return this; }

    public TraversalEngineFactory secondaryColors() { configuration.getSecondaryColors().clear(); return this; }
    public TraversalEngineFactory secondaryColors(int... colors) { Arrays.stream(colors).forEach(c -> configuration.getSecondaryColors().add(c)); return this; }
    public TraversalEngineFactory secondaryColors(Collection<Integer> colors) { configuration.getSecondaryColors().addAll(colors); return this; }

    public TraversalEngineFactory stoppingRule(Class<? extends TraversalStoppingRule<CortexVertex, CortexEdge>> stoppingRule) { configuration.setStoppingRule(stoppingRule); return this; }

    public TraversalEngineFactory graph(DeBruijnGraph graph) { configuration.setGraph(graph); return this; }
    public TraversalEngineFactory rois(DeBruijnGraph rois) { configuration.setRois(rois); return this; }

    public TraversalEngineFactory links() { configuration.getLinks().clear(); return this; }
    public TraversalEngineFactory links(CortexLinks... links) { if (links != null) { Arrays.stream(links).forEach(l -> configuration.getLinks().add(l)); } return this; }
    public TraversalEngineFactory links(Collection<CortexLinks> links) { if (links != null) { configuration.getLinks().addAll(links); } return this; }

    public TraversalEngineFactory references() { configuration.getReferences().clear(); return this; }
    public TraversalEngineFactory references(IndexedReference... lookups) { if (lookups != null) { Arrays.stream(lookups).forEach(r -> configuration.getReferences().add(r)); } return this; }
    public TraversalEngineFactory references(Collection<IndexedReference> lookups) { if (lookups != null) { configuration.getReferences().addAll(lookups); } return this; }

    public TraversalEngineFactory debug() { configuration.setDebugFlag(); return this; }

    public TraversalEngine make() {
        if (configuration.getTraversalColors().size() == 0) {
            throw new CortexJDKException("Traversal color(s) must be specified.");
        }

        for (int c : configuration.getTraversalColors()) {
            if (c >= configuration.getGraph().getNumColors()) {
                throw new CortexJDKException("Traversal colors must be between 0 and " + configuration.getGraph().getNumColors() + " (provided " + c + ")");
            }
        }

        configuration.getJoiningColors().forEach(c -> {
            if (c < 0 || c >= configuration.getGraph().getNumColors()) {
                throw new CortexJDKException("Joining colors must be between 0 and " + configuration.getGraph().getNumColors() + " (provided " + c + ")");
            }
        });

        configuration.getRecruitmentColors().forEach(c -> {
            if (c < 0 || c >= configuration.getGraph().getNumColors()) {
                throw new CortexJDKException("Recruitment colors must be between 0 and " + configuration.getGraph().getNumColors() + " (provided " + c + ")");
            }
        });

        configuration.getSecondaryColors().forEach(c -> {
            if (c < 0 || c >= configuration.getGraph().getNumColors()) {
                throw new CortexJDKException("Secondary colors must be between 0 and " + configuration.getGraph().getNumColors() + " (provided " + c + ")");
            }
        });

        if (configuration.getStoppingRule() == null) { throw new CortexJDKException("Must provide stopping rule for graph traversal"); }

        if (configuration.getGraph() == null) { throw new CortexJDKException("Must provide graph to traverse."); }

        return new TraversalEngine(configuration);
    }
}
