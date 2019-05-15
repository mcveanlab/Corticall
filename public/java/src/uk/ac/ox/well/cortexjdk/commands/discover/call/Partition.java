package uk.ac.ox.well.cortexjdk.commands.discover.call;

import com.google.common.base.Joiner;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.NovelContinuationStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.NovelKmerLimitedContigStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.NovelPartitionStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.REVERSE;

public class Partition extends Module {
    @Argument(fullName = "graph", shortName = "g", doc = "Graph")
    public CortexGraph GRAPH;

    @Argument(fullName = "links", shortName = "l", doc = "Links", required=false)
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName = "roi", shortName = "r", doc = "ROI")
    public CortexGraph ROIS;

    @Argument(fullName = "linkNovels", shortName = "ln", doc = "Link novels")
    public Boolean LINK_NOVELS = false;

    @Argument(fullName = "kmerTable", shortName = "kt", doc = "Kmer table", required = false)
    public File KMER_TABLE;

    @Argument(fullName = "variantTable", shortName = "vt", doc = "Variant table", required = false)
    public File VARIANT_TABLE;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CanonicalKmer, Map<String, String>> kmerMap = new HashMap<>();
        Map<CanonicalKmer, Set<Integer>> kmerIds = new HashMap<>();
        Map<Integer, Map<String, String>> variantMap = new HashMap<>();

        if (KMER_TABLE != null) {
            TableReader ktr = new TableReader(KMER_TABLE, "id", "length", "kmerIndex", "kmer");
            for (Map<String, String> te : ktr) {
                CanonicalKmer ck = new CanonicalKmer(te.get("kmer"));
                kmerMap.put(ck, te);

                if (!kmerIds.containsKey(ck)) {
                    kmerIds.put(ck, new TreeSet<>());
                }

                kmerIds.get(ck).add(Integer.valueOf(te.get("id")));
            }
        }

        if (VARIANT_TABLE != null) {
            TableReader vtr = new TableReader(VARIANT_TABLE);
            for (Map<String, String> te : vtr) {
                int variantId = Integer.valueOf(te.get("index"));
                if (variantId > -1) {
                    variantMap.put(variantId, te);
                }
            }
        }

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColors(getTraversalColor(GRAPH, ROIS))
                .traversalDirection(BOTH)
                .combinationOperator(OR)
                .graph(GRAPH)
                .links(LINKS)
                .rois(ROIS)
                //.stoppingRule(LINK_NOVELS ? NovelPartitionStopper.class : NovelKmerLimitedContigStopper.class)
                .stoppingRule(LINK_NOVELS ? NovelPartitionStopper.class : ContigStopper.class)
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
                //log.info("nr={}", ROIS.findRecord(ck));

                DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = e.dfs(ck);
                List<CortexVertex> w = TraversalUtils.toWalk(g, ck, getTraversalColor(GRAPH, ROIS));

                /*
                int start = -1, stop = w.size();
                for (int i = 0; i < w.size(); i++) {
                    if (used.containsKey(w.get(i).getCanonicalKmer())) {
                        if (start == -1) { start = i - 2000; }
                        stop = i + 2000;
                    }
                }

                if (start < 0) { start = 0; }
                if (stop >= w.size()) { stop = w.size() - 1; }

                w = w.subList(start, stop);
                */

                if (w.size() == 0) {
                    w = new ArrayList<>();
                    w.add(new CortexVertexFactory()
                            .bases(ck.getKmerAsString())
                            .record(GRAPH.findRecord(ck.getKmerAsString()))
                            .make()
                    );
                }

                int numExp = 0, numFound = 0;
                Map<String, String> vm = null;
                if (kmerMap.containsKey(ck)) {
                    int variantId = Integer.valueOf(kmerMap.get(ck).get("id"));
                    if (variantMap.containsKey(variantId)) {
                        vm = variantMap.get(variantId);
                        //log.info("   {}", Joiner.on(", ").withKeyValueSeparator("=").join(variantMap.get(variantId)));

                        String contigWithNewAllele = (variantMap.get(variantId).get("sleft") + (variantMap.get(variantId).get("new").equals(".") ? "" : variantMap.get(variantId).get("new")) + variantMap.get(variantId).get("sright")).toUpperCase();
                        for (int i = 0; i <= contigWithNewAllele.length() - GRAPH.getKmerSize(); i++) {
                            String sk = contigWithNewAllele.substring(i, i + GRAPH.getKmerSize());
                            CanonicalKmer ak = new CanonicalKmer(sk);

                            CortexVertex v = TraversalUtils.findVertex(g, ak);

                            if (used.containsKey(ak)) {
                                numExp++;
                                if (v != null) {
                                    numFound++;
                                }
                            }

                            //log.info("    {}", v);
                        }
                    }
                }

                int numNovelsInSubgraph = countNovels(used, g);
                int subgraphSize = g == null ? 0 : g.vertexSet().size();

                Pair<Integer, Integer> numMarked = markUsedRois(used, w);

                if (numFound != numExp) {
                    /*
                    for (CortexVertex v : w) {
                        log.debug("v={}", v);
                    }
                    */

                    log.info("  * seed={} subgraphSize={} contigSize={} novelsInSubgraph={} novelsNewlyMarked={} novelsPreviouslyMarked={} recovery={}/{}", ck, subgraphSize, w.size(), numNovelsInSubgraph, numMarked.getFirst(), numMarked.getSecond(), numFound, numExp);
                    if (vm == null) {
                        log.info("    {}", "none");
                    } else {
                        log.info("    {}", Joiner.on(", ").withKeyValueSeparator("=").join(vm));
                    }
                }
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
            int numNovels = 0;
            for (int i = 0; i <= partition.length() - GRAPH.getKmerSize(); i++) {
                CanonicalKmer ck = new CanonicalKmer(partition.substring(i, i + GRAPH.getKmerSize()));

                if (used.containsKey(ck)) {
                    numNovels++;
                }
            }

            out.println(">partition" + numPartitions + " len=" + (partition.length() - GRAPH.getKmerSize() + 1) + " numNovels=" + numNovels);
            out.println(partition);

            numPartitions++;
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

