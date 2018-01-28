package uk.ac.ox.well.cortexjdk.commands.call.call;

import org.jgrapht.graph.DirectedWeightedPseudograph;
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
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.NovelContinuationStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

public class CallCandidates extends Module {
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
        TraversalEngine eSimple = configureTraversalEngine(ContigStopper.class);
        TraversalEngine eFurther = configureTraversalEngine(NovelContinuationStopper.class);

        Map<CanonicalKmer, List<CortexVertex>> used = loadRois(ROIS);

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing novel kmers...")
                .message("records processed")
                .maxRecord(used.size())
                .make(log);

        int numContigs = 0;
//        for (CanonicalKmer ck : used.keySet()) {
        for (CanonicalKmer ck : Arrays.asList(new CanonicalKmer("TACTATTGATAATATACAAAATGAAAATGATCATACTATTGATAATA"))) {
            pm.update();

            if (used.get(ck) == null) {
                List<CortexVertex> wFurther = eFurther.walk(ck);
                DirectedWeightedPseudograph<CortexVertex, CortexEdge> gFurther = eFurther.dfs(ck);

                List<CortexVertex> wSimple = eSimple.walk(ck);
                DirectedWeightedPseudograph<CortexVertex, CortexEdge> gSimple = eSimple.dfs(ck);

                int numMarked = markUsedRois(used, wSimple);
                numContigs++;

                log.info("    contig seed={} lenSimple={} dfsSimple={} lenFurther={} dfsFurther={} rois={}", ck, wSimple.size(), gSimple.vertexSet().size(), wFurther.size(), gFurther.vertexSet().size(), numMarked);
            }
        }

        log.info("Assigned {} contigs to {} novel kmers", numContigs, used.size());
    }

    private int markUsedRois(Map<CanonicalKmer, List<CortexVertex>> used, List<CortexVertex> w) {
        int numMarked = 0;
        for (CortexVertex v : w) {
            if (used.containsKey(v.getCanonicalKmer()) && (used.get(v.getCanonicalKmer()) == null || w.size() > used.get(v.getCanonicalKmer()).size())) {
                used.put(v.getCanonicalKmer(), w);
                numMarked++;
            }
        }

        return numMarked;
    }

    private Map<CanonicalKmer, List<CortexVertex>> loadRois(CortexGraph rois) {
        Map<CanonicalKmer, List<CortexVertex>> used = new HashMap<>();
        for (CortexRecord cr : rois) {
            used.put(cr.getCanonicalKmer(), null);
        }

        return used;
    }

    private TraversalEngine configureTraversalEngine(Class stopper) {
        return new TraversalEngineFactory()
                    .traversalColor(getTraversalColor(GRAPH, ROIS))
                    .traversalDirection(BOTH)
                    .combinationOperator(OR)
                    .graph(GRAPH)
                    .links(LINKS)
                    .references(REFERENCES)
                    .rois(ROIS)
                    .stoppingRule(stopper)
                    .make();
    }

    private int getTraversalColor(DeBruijnGraph graph, CortexGraph rois) {
        return graph.getColorForSampleName(rois.getSampleName(0));
    }
}
