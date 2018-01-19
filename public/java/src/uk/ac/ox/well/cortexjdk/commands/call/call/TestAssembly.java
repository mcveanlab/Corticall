package uk.ac.ox.well.cortexjdk.commands.call.call;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
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

    @Argument(fullName = "mc", shortName = "m", doc="McCortex output")
    public FastaSequenceFile MC;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TraversalEngine e = configureTraversalEngine();

        Map<String, Boolean> used = loadRois(MC);

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing novel kmers...")
                .message("records processed")
                .maxRecord(used.size())
                .make(log);

        for (String sk : used.keySet()) {
            List<CortexVertex> g = e.walk(sk);

            String contig = TraversalEngine.toContig(g);

            out.println("seed=" + sk + " len=" + (contig.length() - 46) + " " + contig);

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

    private Map<String, Boolean> loadRois(FastaSequenceFile fa) {
        Map<String, Boolean> used = new HashMap<>();

        ReferenceSequence rseq;
        while ((rseq = fa.nextSequence()) != null) {
            String[] pieces = rseq.getName().split("\\s+");

            for (String piece : pieces) {
                //log.info("{}", piece);

                if (piece.contains("seed=")) {
                    String[] kv = piece.split("=");

                    used.put(kv[1], false);
                }
            }
        }

        return used;
    }
}
