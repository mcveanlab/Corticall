package uk.ac.ox.well.indiana.commands.caller.call;

import com.google.api.client.util.Joiner;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.jgrapht.traverse.TopologicalOrderIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.stoppingconditions.NahrStopper;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;

import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

/**
 * Created by kiran on 05/06/2017.
 */
public class CallNAHRs extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="lookups", shortName="l", doc="Lookups")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(PARENTS);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(GRAPH.getColorForSampleName(CHILD))
                .joiningColors(GRAPH.getColorsForSampleNames(PARENTS))
                .traversalDirection(BOTH)
                .combinationOperator(AND)
                .stopper(NahrStopper.class)
                .graph(GRAPH)
                .rois(ROI)
                .make();

        Map<CortexKmer, Boolean> used = new HashMap<>();
        for (CortexRecord rr : ROI) {
            used.put(rr.getCortexKmer(), false);
        }

        for (CortexKmer rk : used.keySet()) {
            if (!used.get(rk)) {
                log.info("{}", rk);

                DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = e.dfs(rk.getKmerAsString());

                TopologicalOrderIterator<CortexVertex, CortexEdge> toi = new TopologicalOrderIterator<>(g);

                while (toi.hasNext()) {
                    CortexVertex cv = toi.next();

                    for (String parent : LOOKUPS.keySet()) {
                        KmerLookup kl = LOOKUPS.get(parent);

                        Set<Interval> intervals = kl.findKmer(cv.getSk());

                        if (intervals.size() == 1) {
                            log.info(" -- {} {} {} {}", parent, intervals, used.containsKey(cv.getCk()), recordToString(cv.getCr(), childColor, parentColors));
                        } else {
                            log.info(" -- {} {} {} {}", parent, "[many]", used.containsKey(cv.getCk()), recordToString(cv.getCr(), childColor, parentColors));
                        }
                    }

                    if (used.containsKey(cv.getCk())) {
                        used.put(cv.getCk(), true);
                    }
                }
            }
        }
    }

    private String recordToString(CortexRecord cr, int childColor, List<Integer> parentColors) {
        List<Object> pieces = new ArrayList<>();
        pieces.add(cr.getKmerAsString());
        pieces.add(cr.getCoverage(childColor));
        parentColors.forEach(c -> pieces.add(cr.getCoverage(c)));
        pieces.add(cr.getEdgesAsString(childColor));
        parentColors.forEach(c -> pieces.add(cr.getEdgesAsString(c)));

        return Joiner.on(' ').join(pieces);
    }
}
