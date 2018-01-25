package uk.ac.ox.well.cortexjdk.commands.call.call;

import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.NovelContinuationStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;

/**
 * Created by kiran on 03/11/2017.
 */
public class BetterCall extends Module {
    @Argument(fullName = "graph", shortName = "g", doc = "Graph")
    public CortexGraph GRAPH;

    @Argument(fullName = "links", shortName = "l", doc = "Links", required=false)
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName = "references", shortName = "R", doc = "References")
    public HashMap<String, IndexedReference> REFERENCES;

    @Argument(fullName = "roi", shortName = "r", doc = "ROI")
    public CortexGraph ROI;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CanonicalKmer, Boolean> rois = loadRois();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing novel kmers")
                .message("processed")
                .maxRecord(rois.size())
                .make(log);

        for (CanonicalKmer ck : rois.keySet()) {
            if (!rois.get(ck)) {
                List<CortexVertex> lw = longWalk(rois, ck);
                List<CortexVertex> lg = longGraph(rois, ck);

                log.info("{} {} {}", ck, lw.size(), lg.size());
            }

            pm.update();
        }
    }

    private Map<CanonicalKmer, Boolean> loadRois() {
        Map<CanonicalKmer, Boolean> rois = new HashMap<>();
        for (CortexRecord cr : ROI) {
            rois.put(cr.getCanonicalKmer(), false);
        }

        return rois;
    }

    @NotNull
    private List<CortexVertex> longGraph(Map<CanonicalKmer, Boolean> rois, CanonicalKmer ck) {
        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(GRAPH.getColorForSampleName(ROI.getSampleName(0)))
                .joiningColors(GRAPH.getColorsForSampleNames(REFERENCES.keySet()))
                .combinationOperator(OR)
                .stoppingRule(NovelContinuationStopper.class)
                .graph(GRAPH)
                .links(LINKS)
                .references(REFERENCES.values())
                .rois(ROI)
                .make();

        return TraversalEngine.toWalk(e.dfs(ck.getKmerAsString()), ck.getKmerAsString());
    }

    @NotNull
    private List<CortexVertex> longWalk(Map<CanonicalKmer, Boolean> rois, CanonicalKmer ck) {
        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(GRAPH.getColorForSampleName(ROI.getSampleName(0)))
                .joiningColors(GRAPH.getColorsForSampleNames(REFERENCES.keySet()))
                .combinationOperator(OR)
                .stoppingRule(NovelContinuationStopper.class)
                .graph(GRAPH)
                .links(LINKS)
                .references(REFERENCES.values())
                .rois(ROI)
                .make();

        List<CortexVertex> w = e.walk(ck.getKmerAsString());

        boolean extended;
        do {
            extended = false;
            List<List<CortexVertex>> extFwd = new ArrayList<>();

            Set<CortexVertex> nvs = e.getNextVertices(w.get(w.size() - 1).getKmerAsByteKmer());
            for (CortexVertex cv : nvs) {
                List<CortexVertex> wn = e.walk(cv.getKmerAsString(), true);
                wn.add(0, cv);

                boolean hasNovels = false;

                for (CortexVertex v : wn) {
                    if (rois.containsKey(v.getCanonicalKmer()) && !rois.get(v.getCanonicalKmer())) {
                        hasNovels = true;
                        break;
                    }
                }

                if (hasNovels) {
                    extFwd.add(wn);
                }
            }

            if (extFwd.size() == 1) {
                w.addAll(extFwd.get(0));
                extended = true;

                for (CortexVertex v : extFwd.get(0)) {
                    if (rois.containsKey(v.getCanonicalKmer())) {
                        rois.put(v.getCanonicalKmer(), true);
                    }
                }
            }
        } while (extended);

        do {
            extended = false;
            List<List<CortexVertex>> extRev = new ArrayList<>();

            Set<CortexVertex> pvs = e.getPrevVertices(w.get(0).getKmerAsByteKmer());
            for (CortexVertex cv : pvs) {
                List<CortexVertex> wp = e.walk(cv.getKmerAsString(), false);
                wp.add(cv);

                boolean hasNovels = false;

                for (CortexVertex v : wp) {
                    if (rois.containsKey(v.getCanonicalKmer()) && !rois.get(v.getCanonicalKmer())) {
                        hasNovels = true;
                        break;
                    }
                }

                if (hasNovels) {
                    extRev.add(wp);
                }
            }

            if (extRev.size() == 1) {
                w.addAll(0, extRev.get(0));
                extended = true;

                for (CortexVertex v : extRev.get(0)) {
                    if (rois.containsKey(v.getCanonicalKmer())) {
                        rois.put(v.getCanonicalKmer(), true);
                    }
                }
            }
        } while (extended);

        return w;
    }
}
