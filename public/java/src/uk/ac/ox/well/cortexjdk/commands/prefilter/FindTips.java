package uk.ac.ox.well.cortexjdk.commands.prefilter;

import org.jgrapht.DirectedGraph;
import org.jgrapht.Graph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Description;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraphWriter;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.TipBeginningStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.TipEndStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.File;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

@Description(text="Find chains of novel kmers only anchored at one end")
public class FindTips extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName = "links", shortName = "l", doc = "Links", required = false)
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Output
    public File out;

    //@Output(fullName="excluded_out", shortName="xo", doc="Excluded kmers output file")
    //public File tips_out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        Set<Integer> parentColors = new HashSet<>(GRAPH.getColorsForSampleNames(PARENTS));
        //Set<Integer> ignoreColors = new HashSet<>(GRAPH.getColorsForSampleNames(IGNORE));

        log.info("Colors:");
        log.info(" -   child: {}", GRAPH.getColorForSampleName(CHILD));
        log.info(" - parents: {}", GRAPH.getColorsForSampleNames(PARENTS));
        //log.info(" -  ignore: {}", GRAPH.getColorsForSampleNames(IGNORE));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Finding tips")
                .message("records processed")
                .maxRecord(ROI.getNumRecords())
                .make(log);

        Map<CanonicalKmer, Boolean> used = new HashMap<>();
        for (CortexRecord rr : ROI) {
            used.put(rr.getCanonicalKmer(), false);
        }

        Set<CanonicalKmer> tips = new HashSet<>();
        int numTipChains = 0;

        for (CanonicalKmer rr : used.keySet()) {
            if (!used.get(rr)) {
                TraversalEngine e = new TraversalEngineFactory()
                        .traversalDirection(BOTH)
                        .combinationOperator(AND)
                        .traversalColor(childColor)
                        .joiningColors(parentColors)
                        .stoppingRule(ContigStopper.class)
                        .rois(ROI)
                        .graph(GRAPH)
                        .links(LINKS)
                        .make();

                List<CortexVertex> l = e.walk(rr.getKmerAsString());

                boolean novelEnds = used.containsKey(l.get(0).getCanonicalKmer()) || used.containsKey(l.get(l.size() - 1).getCanonicalKmer());

                for (CortexVertex cv : l) {
                    if (used.containsKey(cv.getCanonicalKmer())) {
                        used.put(cv.getCanonicalKmer(), true);

                        if (novelEnds) {
                            tips.add(cv.getCanonicalKmer());
                        }
                    }
                }

                /*
                Graph<CortexVertex, CortexEdge> dfsToParents = null;
                Graph<CortexVertex, CortexEdge> dfsToFree = null;

                for (boolean goForward : Arrays.asList(true, false)) {
                    TraversalEngine eFwd = new TraversalEngineFactory()
                            .traversalDirection(goForward ? FORWARD : REVERSE)
                            .combinationOperator(AND)
                            .traversalColor(childColor)
                            .joiningColors(parentColors)
                            .stoppingRule(TipBeginningStopper.class)
                            .rois(ROI)
                            .graph(GRAPH)
                            .make();

                    TraversalEngine eRev = new TraversalEngineFactory()
                            .traversalDirection(goForward ? REVERSE : FORWARD)
                            .combinationOperator(AND)
                            .traversalColor(childColor)
                            .joiningColors(parentColors)
                            .stoppingRule(TipEndStopper.class)
                            .rois(ROI)
                            .graph(GRAPH)
                            .make();

                    dfsToParents = eFwd.dfs(rr.getKmerAsString());
                    dfsToFree = eRev.dfs(rr.getKmerAsString());

                    if (dfsToParents != null) {
                        for (CortexVertex cv : dfsToParents.vertexSet()) {
                            if (used.containsKey(cv.getCanonicalKmer())) {
                                used.put(cv.getCanonicalKmer(), true);
                            }
                        }
                    }

                    if (dfsToFree != null) {
                        for (CortexVertex cv : dfsToFree.vertexSet()) {
                            if (used.containsKey(cv.getCanonicalKmer())) {
                                used.put(cv.getCanonicalKmer(), true);
                            }
                        }
                    }

                    if (dfsToParents != null && dfsToParents.vertexSet().size() > 0 && dfsToFree != null && dfsToFree.vertexSet().size() > 0) {
                        break;
                    }
                }

                DirectedGraph<CortexVertex, CortexEdge> dfs = null;
                if (dfsToParents != null && dfsToParents.vertexSet().size() > 0 && dfsToFree != null && dfsToFree.vertexSet().size() > 0) {
                    dfs = new DirectedWeightedPseudograph<>(CortexEdge.class);

                    Graphs.addGraph(dfs, dfsToParents);
                    Graphs.addGraph(dfs, dfsToFree);
                }

                if (dfs != null && dfs.vertexSet().size() > 0) {
                    numTipChains++;

                    log.debug("    tip chain {}, seed {}, {} vertices", numTipChains, rr.getKmerAsString(), dfs.vertexSet().size());

                    for (CortexVertex av : dfs.vertexSet()) {
                        tips.add(av.getCanonicalKmer());
                    }
                }
                */
            }

            pm.update();
        }

        log.info("Found {} tip kmer chains ({} kmers total)", numTipChains, tips.size());

        log.info("Writing...");

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(ROI.getHeader());

        //CortexGraphWriter cgt = new CortexGraphWriter(tips_out);
        //cgt.setHeader(ROI.getHeader());

        int numKept = 0, numExcluded = 0;
        for (CortexRecord rr : ROI) {
            if (!tips.contains(rr.getCanonicalKmer())) {
                //cgw.addRecord(rr);
                numKept++;
            } else {
                //cgt.addRecord(rr);
                cgw.addRecord(rr);
                numExcluded++;
            }
        }

        cgw.close();
        //cgt.close();

        log.info("  {}/{} ({}%) kept, {}/{} ({}%) excluded",
                numKept,     ROI.getNumRecords(), 100.0f * (float) numKept / (float) ROI.getNumRecords(),
                numExcluded, ROI.getNumRecords(), 100.0f * (float) numExcluded / (float) ROI.getNumRecords()
        );
    }
}
