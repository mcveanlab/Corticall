package uk.ac.ox.well.cortexjdk.commands.discover.call;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.DestinationStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.NovelPartitionStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.REVERSE;

public class MergePartitions extends Module {
    @Argument(fullName = "partitions", shortName = "p", doc = "Partitions")
    public FastaSequenceFile PARTITIONS;

    @Argument(fullName = "graph", shortName = "g", doc = "Graph")
    public CortexGraph GRAPH;

    @Argument(fullName = "links", shortName = "l", doc = "Links", required=false)
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName = "roi", shortName = "r", doc = "ROI")
    public CortexGraph ROIS;

    @Override
    public void execute() {
        Map<String, Set<ReferenceSequence>> kmerToPartition = new HashMap<>();
        Set<ReferenceSequence> rseqs = new HashSet<>();
        ReferenceSequence rseq;
        while ((rseq = PARTITIONS.nextSequence()) != null) {
            rseqs.add(rseq);

            for (boolean rc : Arrays.asList(false, true)) {
                String seq = rc ? rseq.getBaseString() : SequenceUtils.reverseComplement(rseq.getBaseString());

                String sk0 = seq.substring(0, GRAPH.getKmerSize());
                String sk1 = seq.substring(seq.length() - GRAPH.getKmerSize(), seq.length());

                if (!kmerToPartition.containsKey(sk0)) { kmerToPartition.put(sk0, new HashSet<>()); }
                if (!kmerToPartition.containsKey(sk1)) { kmerToPartition.put(sk1, new HashSet<>()); }

                kmerToPartition.get(sk0).add(rseq);
                kmerToPartition.get(sk1).add(rseq);
            }
        }

        TraversalEngine er = new TraversalEngineFactory()
                .traversalColors(getTraversalColor(GRAPH, ROIS))
                .traversalDirection(REVERSE)
                .combinationOperator(OR)
                .graph(GRAPH)
                .links(LINKS)
                .rois(ROIS)
                .stoppingRule(DestinationStopper.class)
                .make();

        TraversalEngine ef = new TraversalEngineFactory()
                .traversalColors(getTraversalColor(GRAPH, ROIS))
                .traversalDirection(FORWARD)
                .combinationOperator(OR)
                .graph(GRAPH)
                .links(LINKS)
                .rois(ROIS)
                .stoppingRule(DestinationStopper.class)
                .make();

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColors(getTraversalColor(GRAPH, ROIS))
                .traversalDirection(FORWARD)
                .combinationOperator(OR)
                .graph(GRAPH)
                .links(LINKS)
                .stoppingRule(DestinationStopper.class)
                .make();

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = e.dfs("GGAAGTAGTAATTATTGGCAGTTAGTAAGTTCATCATTGTTAAATGC", "TAAAAAAAACCATTTCATCATTAAAAACATATCAAGGAGACATGAAT");

        for (ReferenceSequence rseq1 : rseqs) {
            String sr = rseq1.getBaseString().substring(0, GRAPH.getKmerSize());
            String sf = rseq1.getBaseString().substring(rseq1.length() - GRAPH.getKmerSize(), rseq1.length());

            Set<String> rsinks = new HashSet<>(kmerToPartition.keySet());
            rsinks.remove(sr);
            rsinks.remove(SequenceUtils.reverseComplement(sr));

            Set<String> fsinks = new HashSet<>(kmerToPartition.keySet());
            rsinks.remove(SequenceUtils.reverseComplement(sf));

            DirectedWeightedPseudograph<CortexVertex, CortexEdge> gr = er.dfs(sr, rsinks.toArray(new String[0]));
            List<CortexVertex> wr = TraversalUtils.toWalk(gr, sr, getTraversalColor(GRAPH, ROIS));

            DirectedWeightedPseudograph<CortexVertex, CortexEdge> gf = ef.dfs(sf, fsinks.toArray(new String[0]));
            List<CortexVertex> wf = TraversalUtils.toWalk(gf, sf, getTraversalColor(GRAPH, ROIS));

            log.info("");
        }
    }

    private int getTraversalColor(DeBruijnGraph graph, CortexGraph rois) {
        return graph.getColorForSampleName(rois.getSampleName(0));
    }
}
