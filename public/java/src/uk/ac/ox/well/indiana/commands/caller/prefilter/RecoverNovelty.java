package uk.ac.ox.well.indiana.commands.caller.prefilter;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.stoppingconditions.NovelKmerAggregationStopper;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;

import java.io.File;
import java.util.*;

import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

/**
 * Created by kiran on 26/05/2017.
 */
public class RecoverNovelty extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Output
    public File out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        Set<Integer> parentColors = new HashSet<>(GRAPH.getColorsForSampleNames(PARENTS));

        log.info("Colors:");
        log.info(" -   child: {}", GRAPH.getColorForSampleName(CHILD));
        log.info(" - parents: {}", GRAPH.getColorsForSampleNames(PARENTS));

        TraversalEngine e = new TraversalEngineFactory()
                .traversalDirection(BOTH)
                .combinationOperator(OR)
                .traversalColor(childColor)
                .joiningColors(parentColors)
                .stopper(NovelKmerAggregationStopper.class)
                .graph(GRAPH)
                .rois(ROI)
                .make();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing records")
                .message("records processed")
                .maxRecord(ROI.getNumRecords())
                .make(log);

        Map<CortexKmer, CortexRecord> seen = new TreeMap<>();
        for (CortexRecord rr : ROI) {
            if (!seen.containsKey(rr.getCortexKmer())) {
                log.info("{}", rr);

                DirectedGraph<CortexVertex, CortexEdge> g = e.dfs(rr.getKmerAsString());

                for (CortexVertex cv : g.vertexSet()) {
                    CortexRecord cr = cv.getCr();

                    log.info("  {} {}", g.vertexSet().size(), cr);

                    CortexRecord qr = new CortexRecord(
                            cr.getBinaryKmer(),
                            new int[]{cr.getCoverage(childColor)},
                            new byte[]{cr.getEdges()[childColor]},
                            cr.getKmerSize(),
                            cr.getKmerBits()
                    );

                    seen.put(qr.getCortexKmer(), qr);
                }
            }

            pm.update();
        }

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(ROI.getHeader());

        for (CortexKmer ck : seen.keySet()) {
            cgw.addRecord(seen.get(ck));
        }

        cgw.close();
    }
}
