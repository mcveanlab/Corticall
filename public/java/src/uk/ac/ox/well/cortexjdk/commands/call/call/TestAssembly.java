package uk.ac.ox.well.cortexjdk.commands.call.call;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.containers.ContainerUtils;
import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

public class TestAssembly extends Module {
    @Argument(fullName = "graph", shortName = "g", doc = "Graph")
    public CortexGraph GRAPH;

    @Argument(fullName = "links", shortName = "l", doc = "Links", required=false)
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName = "references", shortName = "R", doc = "References", required=false)
    public ArrayList<IndexedReference> REFERENCES;

    @Argument(fullName = "roi", shortName = "r", doc = "ROI")
    public CortexGraph ROIS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TraversalEngine e = configureTraversalEngine();

        Map<CanonicalKmer, Boolean> used = loadRois(ROIS);

        /*
        Map<CanonicalKmer, Boolean> used = new HashMap<>();
        used.put(new CanonicalKmer("ATAACACTAAAAATTAAATACCAAAAAAAAAAAAAAAAAAAAAAAAT"), false);
        used.put(new CanonicalKmer("ACCCTGAACCCTGAACCCTAAACCCTAAAACCTGAACCCTGAACCCT"), false);
        used.put(new CanonicalKmer("GGTTTAGGGTTTAGGGTTCAGGGTTCAGGGTTCAGGGTTCAGGTTTA"), false);
        */

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing novel kmers...")
                .message("records processed")
                .maxRecord(used.size())
                .make(log);

        for (CanonicalKmer ck : used.keySet()) {
            if (!used.get(ck)) {
                List<CortexVertex> g = e.walk(ck.getKmerAsString());

                out.println(">" + ck);
                out.println(TraversalEngine.toContig(g));

                /*
                List<CortexVertex> gwalk = e.gwalk(ck.getKmerAsString());

                List<CortexVertex> g = e.walk(ck.getKmerAsString(), false);
                CortexVertex v = new CortexVertex(new CortexByteKmer(ck.getKmerAsString()), GRAPH.findRecord(ck));
                g.add(v);

                log.info("{} {} {}", ck, gwalk.size(), g.size());

                Map<String, Integer> cvs = new HashMap<>();
                for (CortexVertex cv : g) {
                    ContainerUtils.increment(cvs, cv.getKmerAsString());
                }

                for (String s : cvs.keySet()) {
                    if (cvs.get(s) > 1) {
                        log.info("  {} {}", s, cvs.get(s));
                    }
                }

                used.put(ck, true);
                */
            }

            pm.update();
        }
    }

    private TraversalEngine configureTraversalEngine() {
        return new TraversalEngineFactory()
                    .traversalColor(getTraversalColor(GRAPH, ROIS))
                    .traversalDirection(BOTH)
                    .combinationOperator(OR)
                    .graph(GRAPH)
                    .links(LINKS)
                    .references(REFERENCES)
                    .rois(ROIS)
                    .stoppingRule(ContigStopper.class)
                    .make();
    }

    private int getTraversalColor(DeBruijnGraph graph, CortexGraph rois) {
        return graph.getColorForSampleName(rois.getSampleName(0));
    }

    private Map<CanonicalKmer, Boolean> loadRois(CortexGraph rois) {
        Map<CanonicalKmer, Boolean> used = new HashMap<>();
        for (CortexRecord cr : rois) {
            used.put(cr.getCanonicalKmer(), false);
        }

        return used;
    }
}
