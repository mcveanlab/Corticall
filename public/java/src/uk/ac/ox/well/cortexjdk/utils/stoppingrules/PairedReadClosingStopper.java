package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by kiran on 02/06/2017.
 */
public class PairedReadClosingStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    private Set<CanonicalKmer> sinks = new HashSet<>();

    @Override
    public boolean hasTraversalSucceeded(TraversalState<CortexVertex> s) {
        if (s.getSinks().size() > 0 && sinks.size() == 0) {
            for (String sink : s.getSinks()) {
                sinks.add(new CanonicalKmer(sink));
            }
        }

        //return s.getSinks().contains(s.getCurrentVertex().getKmerAsString());

        return sinks.contains(s.getCurrentVertex().getCanonicalKmer());
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<CortexVertex> s) {
        return s.getCurrentJunctionDepth() >= 5 || s.getNumAdjacentEdges() == 0 || s.reachedMaxBranchLength();
    }
}
