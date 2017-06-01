package uk.ac.ox.well.indiana.commands.caller.prefilter;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.stoppingconditions.BubbleOpeningStopper;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;
import uk.ac.ox.well.indiana.utils.visualizer.GraphVisualizer;
import uk.ac.ox.well.indiana.utils.visualizer.GraphVisualizationFactory;

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

        GraphVisualizer vc = new GraphVisualizer(9000);

        Set<CortexKmer> errorKmers = new HashSet<>();

        for (CortexRecord rr : ROI) {
            log.info("{}", rr);

            if (!errorKmers.contains(rr.getCortexKmer())) {
                TraversalEngine o = new TraversalEngineFactory()
                        .traversalColor(childColor)
                        .joiningColors(parentColors)
                        .combinationOperator(AND)
                        .traversalDirection(BOTH)
                        .stopper(BubbleOpeningStopper.class)
                        .rois(ROI)
                        .graph(GRAPH)
                        .make();

                Set<Integer> displayColors = new HashSet<>();
                displayColors.add(childColor);
                displayColors.addAll(parentColors);

                DirectedGraph<CortexVertex, CortexEdge> g = o.dfs(rr.getKmerAsString());

                log.info("  {} {}", rr.getCortexKmer(), g.vertexSet().size());

                vc.display(g, displayColors, rr.getKmerAsString());

                for (CortexVertex v : g.vertexSet()) {
                    errorKmers.add(v.getCr().getCortexKmer());
                }

                /*
                if (g1 != null) {
                    TraversalEngine c = new TraversalEngineFactory()
                            .traversalColor(childColor)
                            .joiningColors(parentColors)
                            .combinationOperator(AND)
                            .traversalDirection(FORWARD)
                            .stopper(BubbleClosingStopper.class)
                            .rois(ROI)
                            .graph(GRAPH)
                            .previousTraversal(g1)
                            .make();

                    for (CortexVertex cv : g1.vertexSet()) {
                        if (o.getOutDegree(cv.getSk()) > 1) {
                            DirectedGraph<CortexVertex, CortexEdge> g2 = c.dfs(cv.getSk());

                            if (g2 != null) {
                                DirectedGraph<CortexVertex, CortexEdge> gvis = new DefaultDirectedGraph<>(CortexEdge.class);
                                Graphs.addAllVertices(gvis, g1.vertexSet());
                                Graphs.addAllVertices(gvis, g2.vertexSet());
                                for (CortexEdge e : g1.edgeSet()) {
                                    gvis.addEdge(g1.getEdgeSource(e), g1.getEdgeTarget(e), new CortexEdge(1, 0.0));
                                }
                                for (CortexEdge e : g2.edgeSet()) {
                                    gvis.addEdge(g2.getEdgeSource(e), g2.getEdgeTarget(e), new CortexEdge(1, 0.0));
                                }

                                VisualCortex vc = new VisualCortexFactory()
                                        .subgraph(gvis)
                                        .port(9000)
                                        .logger(log)
                                        .make();

                                DepthFirstIterator<CortexVertex, CortexEdge> dfi1 = new DepthFirstIterator<>(g1, cv);
                                DepthFirstIterator<CortexVertex, CortexEdge> dfi2 = new DepthFirstIterator<>(g2, cv);

                                Set<GraphPath<CortexVertex, CortexEdge>> ps1 = new HashSet<>();

                                Set<CortexVertex> ends = new HashSet<>();
                                while (dfi1.hasNext()) {
                                    CortexVertex ev = dfi1.next();

                                    if (o.getInDegree(ev.getSk()) > 1) {
                                        ends.add(ev);

                                        DijkstraShortestPath<CortexVertex, CortexEdge> dj = new DijkstraShortestPath<>(g1, cv, ev);

                                        ps1.add(dj.getPath());
                                    }
                                }

                                Set<GraphPath<CortexVertex, CortexEdge>> ps2 = new HashSet<>();

                                while (dfi2.hasNext()) {
                                    CortexVertex ev = dfi2.next();

                                    if (ends.contains(ev)) {
                                        DijkstraShortestPath<CortexVertex, CortexEdge> dj = new DijkstraShortestPath<>(g1, cv, ev);

                                        if (!ps1.contains(dj.getPath())) {
                                            ps2.add(dj.getPath());
                                        }
                                    }
                                }

                                for (GraphPath<CortexVertex, CortexEdge> p1 : ps1) {
                                    float covMean1 = 0.0f;
                                    for (CortexVertex c1 : p1.getVertexList()) {
                                        covMean1 += c1.getCr().getCoverage(childColor);
                                    }
                                    covMean1 = covMean1 / p1.getLength();

                                    for (GraphPath<CortexVertex, CortexEdge> p2 : ps2) {
                                        float covMean2 = 0.0f;
                                        for (CortexVertex c2 : p2.getVertexList()) {
                                            covMean2 += c2.getCr().getCoverage(childColor);
                                        }
                                        covMean2 = covMean2 / p2.getLength();

                                        log.info(" - {} {} {} {}", covMean1, covMean2, p1.getLength(), p2.getLength());
                                    }
                                }
                            }
                        }
                    }
                }
                */
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
