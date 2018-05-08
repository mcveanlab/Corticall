package uk.ac.ox.well.cortexjdk.utils.stoppingrules;

import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalState;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by kiran on 10/05/2017.
 */
public class NovelKmerLimitedContigStopper extends AbstractTraversalStoppingRule<CortexVertex, CortexEdge> {
    private boolean foundNovelKmers = false;
    private int distanceFromSeed = 0;

    private Set<CanonicalKmer> rois = null;

    @Override
    public boolean hasTraversalSucceeded(TraversalState<CortexVertex> s) {
        distanceFromSeed++;

        if (s.getRois() != null) {
            if (rois == null) {
                rois = new HashSet<>();
                for (CortexRecord cr : s.getRois()) {
                    rois.add(cr.getCanonicalKmer());
                }
            }
        } else {
            throw new CortexJDKException("This stopper requires a list of novel kmers be provided.");
        }

        if (rois != null && rois.contains(s.getCurrentVertex().getCanonicalKmer())) {
            foundNovelKmers = true;
            distanceFromSeed = 0;
        }

        boolean stopNow = distanceFromSeed > 2000 || s.getNumAdjacentEdges() != 1 || s.reachedMaxBranchLength();

        return foundNovelKmers && stopNow;
    }

    @Override
    public boolean hasTraversalFailed(TraversalState<CortexVertex> s) {
        return false;
    }
}
