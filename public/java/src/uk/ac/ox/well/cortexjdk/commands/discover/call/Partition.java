package uk.ac.ox.well.cortexjdk.commands.discover.call;

import org.apache.commons.math3.util.Pair;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.NovelContinuationStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

public class Partition extends Module {
    @Argument(fullName = "graph", shortName = "g", doc = "Graph")
    public CortexGraph GRAPH;

    @Argument(fullName = "links", shortName = "l", doc = "Links", required=false)
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName = "roi", shortName = "r", doc = "ROI")
    public CortexGraph ROIS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TraversalEngine e = new TraversalEngineFactory()
                .traversalColors(getTraversalColor(GRAPH, ROIS))
                .traversalDirection(BOTH)
                .combinationOperator(OR)
                .graph(GRAPH)
                .links(LINKS)
                .stoppingRule(ContigStopper.class)
                .make();

        log.info("Using stopper {}", e.getConfiguration().getStoppingRule().getSimpleName());

        Map<CanonicalKmer, List<CortexVertex>> used = loadRois(ROIS);

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing novel kmers...")
                .message("records processed")
                .maxRecord(used.size())
                .make(log);

        for (CanonicalKmer ck : used.keySet()) {
            pm.update();

            if (used.get(ck) == null) {
                DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = e.dfs(ck);
                List<CortexVertex> w = TraversalUtils.toWalk(g, ck, getTraversalColor(GRAPH, ROIS));

                if (w.size() == 0) {
                    w = new ArrayList<>();
                    w.add(new CortexVertexFactory()
                            .bases(ck.getKmerAsString())
                            .record(GRAPH.findRecord(ck.getKmerAsString()))
                            .make()
                    );
                }

                int numNovelsInSubgraph = countNovels(used, g);
                int subgraphSize = g == null ? 0 : g.vertexSet().size();

                Pair<Integer, Integer> numMarked = markUsedRois(used, w);

                log.info("  * seed={} subgraphSize={} contigSize={} novelsInSubgraph={} novelsNewlyMarked={} novelsPreviouslyMarked={}", ck, subgraphSize, w.size(), numNovelsInSubgraph, numMarked.getFirst(), numMarked.getSecond());
            }
        }

        Set<String> contigs = new TreeSet<>();

        int numNovelKmersAssigned = 0;
        for (CanonicalKmer ck : used.keySet()) {
            if (used.get(ck) != null) {
                String fw = TraversalUtils.toContig(used.get(ck));
                String rc = SequenceUtils.reverseComplement(fw);

                if (!contigs.contains(fw) && !contigs.contains(rc)) {
                    contigs.add(fw);
                }

                numNovelKmersAssigned++;
            }
        }

        int numPartitions = 0;
        for (String partition : contigs) {
            out.println(">partition" + numPartitions);
            out.println(partition);
        }

        log.info("Assigned {}/{} novel kmers to {} contigs", numNovelKmersAssigned, used.size(), contigs.size());
    }

    private int countNovels(Map<CanonicalKmer, List<CortexVertex>> used, DirectedWeightedPseudograph<CortexVertex, CortexEdge> g) {
        int numNovels = 0;

        if (g != null) {
            for (CortexVertex v : g.vertexSet()) {
                if (used.containsKey(v.getCanonicalKmer())) {
                    numNovels++;
                }
            }
        }

        return numNovels;
    }

    private Pair<Integer, Integer> markUsedRois(Map<CanonicalKmer, List<CortexVertex>> used, List<CortexVertex> w) {
        int numNewlyMarked = 0, numAlreadyMarked = 0;
        for (CortexVertex v : w) {
            if (used.containsKey(v.getCanonicalKmer()) && used.get(v.getCanonicalKmer()) != null) {
                numAlreadyMarked++;
            }
        }

        for (CortexVertex v : w) {
            if (used.containsKey(v.getCanonicalKmer())) {
                if ((used.get(v.getCanonicalKmer()) == null || w.size() > used.get(v.getCanonicalKmer()).size())) {
                    used.put(v.getCanonicalKmer(), w);
                    numNewlyMarked++;
                }
            }
        }

        return new Pair<>(numNewlyMarked, numAlreadyMarked);
    }

    private Map<CanonicalKmer, List<CortexVertex>> loadRois(CortexGraph rois) {
        Map<CanonicalKmer, List<CortexVertex>> used = new TreeMap<>();
        for (CortexRecord cr : rois) {
            used.put(cr.getCanonicalKmer(), null);
        }

        return used;
    }

    private int getTraversalColor(DeBruijnGraph graph, CortexGraph rois) {
        return graph.getColorForSampleName(rois.getSampleName(0));
    }
}

