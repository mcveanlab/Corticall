package uk.ac.ox.well.indiana.commands.caller.prefilter;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.stoppingconditions.ExplorationStopper;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;
import uk.ac.ox.well.indiana.utils.visualizer.GraphVisualizer;

import java.util.*;

import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

public class RemoveSequencingErrors extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROIs")
    public CortexGraph ROI;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        Set<Integer> parentColors = new HashSet<>(GRAPH.getColorsForSampleNames(PARENTS));

        log.info("Colors:");
        log.info(" -   child: {}", GRAPH.getColorForSampleName(CHILD));
        log.info(" - parents: {}", GRAPH.getColorsForSampleNames(PARENTS));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Finding sequencing errors")
                .message("records processed")
                .maxRecord(ROI.getNumRecords())
                .make(log);

        Set<CortexKmer> roiKmers = new HashSet<>();
        for (CortexRecord rr : ROI) {
            roiKmers.add(rr.getCortexKmer());
        }

        GraphVisualizer vc = new GraphVisualizer(9000);

        Set<CortexKmer> errorKmers = new HashSet<>();

        for (CortexRecord rr : ROI) {
            //log.info("{}", rr);

            if (!errorKmers.contains(rr.getCortexKmer())) {
                TraversalEngine o = new TraversalEngineFactory()
                        .traversalColor(childColor)
                        .joiningColors(parentColors)
                        .traversalDirection(BOTH)
                        .combinationOperator(AND)
                        .connectAllNeighbors(true)
                        //.stopper(BubbleOpeningStopper.class)
                        .stopper(ExplorationStopper.class)
                        .rois(ROI)
                        .graph(GRAPH)
                        .make();

                Set<Integer> displayColors = new HashSet<>();
                displayColors.add(childColor);
                displayColors.addAll(parentColors);

                DirectedGraph<CortexVertex, CortexEdge> g = o.dfs(rr.getKmerAsString());

                int numNovels = 0;
                for (CortexVertex cv : g.vertexSet()) {
                    if (roiKmers.contains(cv.getCk())) {
                        numNovels++;
                    }
                }

                log.info("  rr: {}, vertices: {}, novels: {}", rr, g.vertexSet().size(), numNovels);

                vc.display(g, rr.getKmerAsString());

                for (CortexVertex v : g.vertexSet()) {
                    if (roiKmers.contains(v.getCk())) {
                        errorKmers.add(v.getCr().getCortexKmer());
                    }
                }
            }

            pm.update();
        }

        log.info("Found {} sequencing errors", errorKmers.size());

        /*
        log.info("Writing...");

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(ROI.getHeader());

        CortexGraphWriter cgo = new CortexGraphWriter(shared_out);
        cgo.setHeader(ROI.getHeader());

        int numKept = 0, numExcluded = 0;
        for (CortexRecord rr : ROI) {
            if (!sharedKmers.contains(rr.getCortexKmer())) {
                cgw.addRecord(rr);
                numKept++;
            } else {
                cgo.addRecord(rr);
                numExcluded++;
            }
        }

        cgw.close();
        cgo.close();

        log.info("  {}/{} ({}%) kept, {}/{} ({}%) excluded",
                numKept,     ROI.getNumRecords(), 100.0f * (float) numKept / (float) ROI.getNumRecords(),
                numExcluded, ROI.getNumRecords(), 100.0f * (float) numExcluded / (float) ROI.getNumRecords()
        );
        */
    }
}
