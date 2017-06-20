package uk.ac.ox.well.indiana.commands.caller.call;

import htsjdk.samtools.util.Interval;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.jgrapht.graph.EdgeReversedGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksMap;
import uk.ac.ox.well.indiana.utils.stoppingconditions.ContigStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.NahrStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.TraversalStopper;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;

import java.util.*;

import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

/**
 * Created by kiran on 20/06/2017.
 */
public class CallNahrEvents extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Links")
    public CortexLinksMap LINKS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="refs", shortName="R", doc="References")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(PARENTS);
        List<Integer> recruitColors = GRAPH.getColorsForSampleNames(new ArrayList<>(LOOKUPS.keySet()));

        logColorAssignments(childColor, parentColors, recruitColors);

        Map<CortexKmer, Boolean> usedRois = loadRois();

        TraversalEngine ce = initializeTraversalEngine(childColor, parentColors, recruitColors, ContigStopper.class);
        //TraversalEngine ne = initializeTraversalEngine(childColor, parentColors, recruitColors, NahrStopper.class);

        for (CortexRecord rr : ROI) {
            if (!usedRois.get(rr.getCortexKmer())) {
                DirectedWeightedPseudograph<CortexVertex, CortexEdge> cg = ce.dfs(rr.getKmerAsString());
                CortexVertex rv = new CortexVertex(rr.getKmerAsString(), GRAPH.findRecord(rr));

                for (String key : LOOKUPS.keySet()) {
                    KmerLookup kl = LOOKUPS.get(key);

                    DepthFirstIterator<CortexVertex, CortexEdge> dfsf = new DepthFirstIterator<>(cg, rv);
                    Map<String, Integer> fContigCount = getContigCounts(kl, dfsf);

                    DepthFirstIterator<CortexVertex, CortexEdge> dfsr = new DepthFirstIterator<>(new EdgeReversedGraph<>(cg), rv);
                    Map<String, Integer> rContigCount = getContigCounts(kl, dfsr);

                    log.info("fContigCount: {} {}", key, fContigCount);
                    log.info("rContigCount: {} {}", key, rContigCount);
                }

                for (CortexVertex cv : cg.vertexSet()) {
                    usedRois.put(cv.getCk(), true);
                }
            }
        }
    }

    private Map<String, Integer> getContigCounts(KmerLookup kl, DepthFirstIterator<CortexVertex, CortexEdge> dfs) {
        Map<String, Integer> contigCounts = new HashMap<>();
        while (dfs.hasNext()) {
            CortexVertex cv = dfs.next();
            Set<Interval> loci = kl.findKmer(cv.getSk());
            for (Interval locus : loci) {
                ContainerUtils.increment(contigCounts, locus.getContig());
            }
        }

        return contigCounts;
    }

    private TraversalEngine initializeTraversalEngine(int childColor, List<Integer> parentColors, List<Integer> recruitColors, Class<? extends TraversalStopper<CortexVertex, CortexEdge>> stoppingRule) {
        TraversalEngine e = new TraversalEngineFactory()
            .combinationOperator(AND)
            .traversalDirection(BOTH)
            .traversalColor(childColor)
            .joiningColors(parentColors)
            .recruitmentColors(recruitColors)
            .rois(ROI)
            .connectAllNeighbors(true)
            .stopper(stoppingRule)
            .graph(GRAPH)
            .make();

        return e;
    }

    private Map<CortexKmer, Boolean> loadRois() {
        Map<CortexKmer, Boolean> usedRois = new HashMap<>();
        for (CortexRecord rr : ROI) {
            usedRois.put(rr.getCortexKmer(), false);
        }

        return usedRois;
    }

    private void logColorAssignments(int childColor, List<Integer> parentColors, List<Integer> recruitColors) {
        log.info("Colors:");
        log.info("    child: {} {}", childColor, CHILD);
        for (int i = 0; i < parentColors.size(); i++) {
            log.info("   parent: {} {}", parentColors.get(i), GRAPH.getSampleName(parentColors.get(i)));
        }
        for (int i = 0; i < recruitColors.size(); i++) {
            log.info("  recruit: {} {}", recruitColors.get(i), GRAPH.getSampleName(recruitColors.get(i)));
        }
    }
}
