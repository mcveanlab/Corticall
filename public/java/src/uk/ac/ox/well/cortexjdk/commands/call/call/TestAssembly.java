package uk.ac.ox.well.cortexjdk.commands.call.call;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

public class TestAssembly extends Module {
    @Argument(fullName = "graph", shortName = "g", doc = "Graph")
    public CortexGraph GRAPH;

    @Argument(fullName = "links", shortName = "l", doc = "Links", required=false)
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName = "mc", shortName = "m", doc="McCortex output")
    public FastaSequenceFile MC;

    @Argument(fullName = "useDfs", shortName = "u", doc="DFS")
    public Boolean USE_DFS = false;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TraversalEngine e = configureTraversalEngine();

//        Map<String, String> used = loadRois(MC);

        Map<String, String> usedFull = loadRois(MC);
        Map<String, String> used = new HashMap<>();
        used.put("ACTTCTACTAACTTAGGTCTTACATTAACTAACTCAGGTCTTACTTT", usedFull.get("ACTTCTACTAACTTAGGTCTTACATTAACTAACTCAGGTCTTACTTT"));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing novel kmers...")
                .message("records processed")
                .maxRecord(used.size())
                .make(log);

        for (String sk : used.keySet()) {
            List<CortexVertex> gw = e.walk(sk);
            List<CortexVertex> gd = e.gwalk(sk);

            String cd = TraversalEngine.toContig(gd);
            String cw = TraversalEngine.toContig(gw);

            String contig = USE_DFS ? cd : cw;

            out.println("seed=" + sk + " len=" + (contig.length() - (GRAPH.getKmerSize() - 1)) + " " + contig);

            log.info("{} cjd={} cjw={} mc={}", sk, cd.length(), cw.length(), used.get(sk).length());

            pm.update();
        }
    }

    private TraversalEngine configureTraversalEngine() {
        return new TraversalEngineFactory()
                    .traversalColor(0)
                    .traversalDirection(BOTH)
                    .combinationOperator(OR)
                    .graph(GRAPH)
                    .links(LINKS)
                    .stoppingRule(ContigStopper.class)
                    .make();
    }

    private Map<String, String> loadRois(FastaSequenceFile fa) {
        Map<String, String> used = new LinkedHashMap<>();

        ReferenceSequence rseq;
        while ((rseq = fa.nextSequence()) != null) {
            String[] pieces = rseq.getName().split("\\s+");

            for (String piece : pieces) {
                if (piece.contains("seed=")) {
                    String[] kv = piece.split("=");

                    used.put(kv[1], rseq.getBaseString());
                }
            }
        }

        return used;
    }
}
