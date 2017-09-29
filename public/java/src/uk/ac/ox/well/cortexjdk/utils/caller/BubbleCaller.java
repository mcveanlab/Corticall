package uk.ac.ox.well.cortexjdk.utils.caller;

import htsjdk.samtools.util.Interval;
import org.jgrapht.GraphPath;
import org.jgrapht.alg.shortestpath.DijkstraShortestPath;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.jgrapht.graph.EdgeReversedGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.BubbleClosingStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.BubbleOpeningStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;

/**
 * Created by kiran on 30/08/2017.
 */
public class BubbleCaller {
    private BubbleCallerConfiguration bc;

    private TraversalEngine eOpen;
    private Map<Integer, TraversalEngine> eCloses = new HashMap<>();

    private Set<CortexKmer> rois = new HashSet<>();

    public BubbleCaller(BubbleCallerConfiguration bc) {
        this.bc = bc;

        eOpen = new TraversalEngineFactory()
                .traversalColor(bc.getAlternateColor())
                .joiningColors(bc.getReferenceColors())
                .traversalDirection(BOTH)
                .combinationOperator(AND)
                .stoppingRule(BubbleOpeningStopper.class)
                .graph(bc.getGraph())
                .links(bc.getLinks())
                .rois(bc.getRois())
                .make();

        for (int pc : bc.getReferenceColors()) {
            TraversalEngine eClose = new TraversalEngineFactory()
                    .traversalColor(pc)
                    .joiningColors(bc.getReferenceColors())
                    .traversalDirection(FORWARD)
                    .combinationOperator(OR)
                    .stoppingRule(BubbleClosingStopper.class)
                    .graph(bc.getGraph())
                    .links(bc.getLinks())
                    .rois(bc.getRois())
                    .make();

            eCloses.put(pc, eClose);
        }

        for (CortexRecord rr : bc.getRois()) {
            rois.add(rr.getCortexKmer());
        }
    }

    public Set<Bubble> call(String seed) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gc = eOpen.dfs(seed);

        Set<Bubble> bubbles = new HashSet<>();
        if (gc != null) {
            DepthFirstIterator<CortexVertex, CortexEdge> dFwd = null;
            DepthFirstIterator<CortexVertex, CortexEdge> dRev = null;

            for (CortexVertex cv : gc.vertexSet()) {
                if (cv.getSk().equals(seed)) {
                    dFwd = new DepthFirstIterator<>(gc, cv);
                    dRev = new DepthFirstIterator<>(new EdgeReversedGraph<>(gc), cv);
                }
            }

            Set<CortexVertex> sources = getCandidates(dRev);
            Set<CortexVertex> sinks = getCandidates(dFwd);

            for (CortexVertex so : sources) {
                for (CortexVertex si : sinks) {
                    DirectedWeightedPseudograph<CortexVertex, CortexEdge> s = new DirectedWeightedPseudograph<>(CortexEdge.class);
                    s.addVertex(si);

                    for (int pc : bc.getReferenceColors()) {
                        Set<Interval> soIntervals = bc.getReferences().get(bc.getGraph().getSampleName(pc)).find(so.getSk());
                        Set<Interval> siIntervals = bc.getReferences().get(bc.getGraph().getSampleName(pc)).find(si.getSk());

                        if (so.getCr().getCoverage(pc) > 0 && si.getCr().getCoverage(pc) > 0 && soIntervals.size() == 1 && siIntervals.size() == 1) {
                            Interval soInterval = soIntervals.iterator().next();
                            Interval siInterval = siIntervals.iterator().next();

                            if (soInterval.getContig().equals(siInterval.getContig())) {
                                eCloses.get(pc).getConfiguration().setPreviousTraversal(s);

                                DirectedWeightedPseudograph<CortexVertex, CortexEdge> gp = eCloses.get(pc).dfs(so.getSk());

                                if (gp != null) {
                                    GraphPath<CortexVertex, CortexEdge> dgc = DijkstraShortestPath.findPathBetween(gc, so, si);
                                    GraphPath<CortexVertex, CortexEdge> dgp = DijkstraShortestPath.findPathBetween(gp, so, si);

                                    Set<CortexKmer> novelKmersInBubble = new HashSet<>();
                                    for (CortexVertex cv : dgc.getVertexList()) {
                                        if (rois.contains(cv.getCk())) {
                                            novelKmersInBubble.add(cv.getCk());
                                        }
                                    }

                                    Bubble b = new Bubble(dgp, dgc, bc.getReferences().get(bc.getGraph().getSampleName(pc)), novelKmersInBubble);

                                    bubbles.add(b);
                                }
                            }
                        }
                    }
                }
            }
        }

        return bubbles;
    }

    private Set<CortexVertex> getCandidates(DepthFirstIterator<CortexVertex, CortexEdge> d) {
        Set<CortexVertex> candidates = new HashSet<>();
        if (d != null) {
            while (d.hasNext()) {
                CortexVertex cv = d.next();
                Set<Interval> allIntervals = new HashSet<>();
                for (KmerLookup kl : bc.getReferences().values()) {
                    Set<Interval> intervals = kl.find(cv.getSk());
                    if (intervals != null) {
                        allIntervals.addAll(intervals);
                    }
                }

                if (allIntervals.size() == 1) {
                    candidates.add(cv);
                }
            }
        }
        return candidates;
    }
}
