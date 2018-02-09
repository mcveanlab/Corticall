package uk.ac.ox.well.cortexjdk.utils.alignment.graph;

import htsjdk.samtools.fastq.FastqRecord;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.GraphPath;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.PairedReadClosingStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

public class ReadToGraphAligner {
    private CortexGraph graph;
    private List<CortexLinks> links;
    private TraversalEngine[] engines;
    private List<Integer> colors = new ArrayList<>();
    private Map<Integer, Map<Pair<CanonicalKmer, CanonicalKmer>, Set<Integer>>> cachedDists = new HashMap<>();

    private DirectedWeightedPseudograph<CortexVertex, CortexEdge> sg = new DirectedWeightedPseudograph<>(CortexEdge.class);

    public ReadToGraphAligner(CortexGraph g, ArrayList<CortexLinks> l, int... colors) {
        this.graph = g;
        this.links = l;

        this.engines = new TraversalEngine[g.getNumColors()];
        for (int c : colors) {
            this.colors.add(c);

            engines[c] = new TraversalEngineFactory()
                    .traversalColor(c)
                    .traversalDirection(BOTH)
                    .combinationOperator(OR)
                    .maxBranchLength(5000)
                    //.connectAllNeighbors(true)
                    .secondaryColors(colors)
                    .stoppingRule(PairedReadClosingStopper.class)
                    .graph(graph)
                    .links(links)
                    .make();

            cachedDists.put(c, new HashMap<>());
        }
    }

    public ReadToGraphAligner(CortexGraph g, ArrayList<CortexLinks> l, Collection<Integer> colors) {
        this.graph = g;
        this.links = l;

        this.engines = new TraversalEngine[g.getNumColors()];
        for (int c : colors) {
            this.colors.add(c);

            engines[c] = new TraversalEngineFactory()
                    .traversalColor(c)
                    .traversalDirection(BOTH)
                    .combinationOperator(OR)
                    .maxBranchLength(5000)
                    //.connectAllNeighbors(true)
                    .secondaryColors(colors)
                    .stoppingRule(PairedReadClosingStopper.class)
                    .graph(graph)
                    .links(links)
                    .make();

            cachedDists.put(c, new HashMap<>());
        }
    }

    public SingleEndAlignmentInfo align(FastqRecord fq) {
        String[] sks = new String[fq.length() - graph.getKmerSize() + 1];
        CortexRecord[] crs = new CortexRecord[fq.length() - graph.getKmerSize() + 1];

        for (int i = 0; i <= fq.length() - graph.getKmerSize(); i++) {
            byte[] bk = fq.getReadString().substring(i, i + graph.getKmerSize()).getBytes();
            CortexRecord cr = graph.findRecord(bk);

            sks[i] = new String(bk);
            crs[i] = cr;
        }

        return new SingleEndAlignmentInfo(fq, sks, crs);
    }

    public PairedEndAlignmentInfo align2(FastqRecord fq1, FastqRecord fq2) {
        Set<String> cks1 = new HashSet<>();
        Set<String> cks2 = new HashSet<>();

        for (int i = 0; i <= fq1.length() - graph.getKmerSize(); i++) {
            cks1.add(fq1.getReadString().substring(i, i + graph.getKmerSize()));
        }

        for (int i = 0; i <= fq2.length() - graph.getKmerSize(); i++) {
            cks2.add(fq2.getReadString().substring(i, i + graph.getKmerSize()));
        }

        for (int c : colors) {
        }

        return null;
    }

    public PairedEndAlignmentInfo align(FastqRecord fq1, FastqRecord fq2) {
        SingleEndAlignmentInfo sai1 = align(fq1);
        SingleEndAlignmentInfo sai2 = align(fq2);

        Map<CanonicalKmer, Integer> distFrom5pEnd = new HashMap<>();
        for (int i = 0; i < sai1.getRecords().length; i++) {
            if (sai1.getRecords()[i] != null) {
                distFrom5pEnd.put(sai1.getRecords()[i].getCanonicalKmer(), i);
            }
        }

        for (int i = sai2.getRecords().length - 1; i >= 0; i--) {
            if (sai2.getRecords()[i] != null) {
                distFrom5pEnd.put(sai2.getRecords()[i].getCanonicalKmer(), i);
            }
        }

        Map<Integer, Set<Integer>> allDists = new HashMap<>();

        Map<Integer, Map<Pair<CanonicalKmer, CanonicalKmer>, List<GraphPath<CortexVertex, CortexEdge>>>> cachedPaths = new HashMap<>();
        for (int q = 0; q < colors.size(); q++) {
            cachedPaths.put(colors.get(q), new HashMap<>());
        }

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = null;

        for (int q = 0; q < colors.size(); q++) {
            int c = colors.get(q);
            allDists.put(c, new HashSet<>());

            List<String> k1 = sai1.getAvailableKmers(c);
            List<String> k2 = reverseComplement(sai2.getAvailableKmers(c));

            if (k1.size() > 0 && k2.size() > 0) {
                String sk1 = k1.get(k1.size() - 1);
                String sk2 = k2.get(0);

                CanonicalKmer ck1 = new CanonicalKmer(sk1);
                CanonicalKmer ck2 = new CanonicalKmer(sk2);

                Pair<CanonicalKmer, CanonicalKmer> key = ck1.compareTo(ck2) <= 0 ? new Pair<>(ck1, ck2) : new Pair<>(ck2, ck1);

                if (!cachedDists.get(c).containsKey(key)) {
                    if (g == null) {
                        TraversalEngine ef = new TraversalEngineFactory().configuration(engines[c].getConfiguration()).sink(sk2).make();
                        g = ef.dfs(sk1);

                        if (g == null) {
                            TraversalEngine eb = new TraversalEngineFactory().configuration(engines[c].getConfiguration()).sink(sk1).make();
                            g = eb.dfs(sk2);
                        }
                    } else {
                        Set<String> newOuts = new HashSet<>();
                        Set<String> newIns = new HashSet<>();

                        for (CortexVertex v : g.vertexSet()) {
                            boolean isNotFlipped = v.getKmerAsString().equals(v.getCanonicalKmer().getKmerAsString());

                            Collection<String> remainingInEdges  = isNotFlipped ? v.getCortexRecord().getInEdgesAsStrings(c, false)  : v.getCortexRecord().getOutEdgesAsStrings(c, true);
                            Collection<String> remainingOutEdges = isNotFlipped ? v.getCortexRecord().getOutEdgesAsStrings(c, false) : v.getCortexRecord().getInEdgesAsStrings(c, true);
                            for (int r = 0; r < q; r++) {
                                Collection<String> prevInEdges   = isNotFlipped ? v.getCortexRecord().getInEdgesAsStrings(colors.get(r), false)  : v.getCortexRecord().getOutEdgesAsStrings(colors.get(r), true);
                                Collection<String> prevOutEdges  = isNotFlipped ? v.getCortexRecord().getOutEdgesAsStrings(colors.get(r), false) : v.getCortexRecord().getInEdgesAsStrings(colors.get(r), true);

                                remainingInEdges.removeAll(prevInEdges);
                                remainingOutEdges.removeAll(prevOutEdges);
                            }

                            for (String roe : remainingOutEdges) {
                                newOuts.add(v.getKmerAsString().substring(1, v.getKmerAsString().length()) + roe);
                            }

                            for (String rie : remainingInEdges) {
                                newIns.add(rie + v.getKmerAsString().substring(0, v.getKmerAsString().length() - 1));
                            }
                        }

                        for (String newOut : newOuts) {
                            TraversalEngine ef = new TraversalEngineFactory().configuration(engines[c].getConfiguration()).sink(newIns).make();
                            DirectedWeightedPseudograph<CortexVertex, CortexEdge> g1 = ef.dfs(newOut);

                            if (g1 != null) {
                                Graphs.addGraph(g, g1);
                            }
                        }

                        for (String newIn : newIns) {
                            TraversalEngine ef = new TraversalEngineFactory().configuration(engines[c].getConfiguration()).sink(newOuts).make();
                            DirectedWeightedPseudograph<CortexVertex, CortexEdge> g1 = ef.dfs(newIn);

                            if (g1 != null) {
                                Graphs.addGraph(g, g1);
                            }
                        }
                    }

                    if (g != null) {
                        CortexVertex vSource = TraversalEngine.findVertex(g, ck1);
                        CortexVertex vSink = TraversalEngine.findVertex(g, ck2);

                        if (vSource != null && vSink != null) {
                            List<GraphPath<CortexVertex, CortexEdge>> ps = cachedPaths.get(c).containsKey(key) ? cachedPaths.get(c).get(key) : new PathFinder(g, c).getPaths(vSource, vSink);

                            for (int qc : colors) {
                                boolean isIdentical = true;
                                if (qc != c) {
                                    for (CortexVertex v : g.vertexSet()) {
                                        if (v.getCortexRecord().getEdges()[qc] != v.getCortexRecord().getEdges()[c]) {
                                            isIdentical = false;
                                            break;
                                        }
                                    }
                                }

                                if (isIdentical) {
                                    cachedPaths.get(c).put(key, ps);
                                }
                            }

                            Set<Integer> d = new TreeSet<>();
                            for (GraphPath<CortexVertex, CortexEdge> p : ps) {
                                d.add(p.getLength());
                            }

                            cachedDists.get(c).put(key, d);

                            for (int i = 0; i < sai1.getRecords().length; i++) {
                                for (int j = sai2.getRecords().length - 1; j >= 0; j--) {
                                    if (sai1.getRecords()[i] != null && sai2.getRecords()[j] != null) {
                                        CanonicalKmer tck1 = sai1.getRecords()[i].getCanonicalKmer();
                                        CanonicalKmer tck2 = sai2.getRecords()[j].getCanonicalKmer();

                                        if (!(ck1.equals(tck1) && ck2.equals(tck2))) {
                                            Pair<CanonicalKmer, CanonicalKmer> newkey = tck1.compareTo(tck2) <= 0 ? new Pair<>(tck1, tck2) : new Pair<>(tck2, tck1);

                                            if (!cachedDists.get(c).containsKey(newkey)) {
                                                cachedDists.get(c).put(newkey, new TreeSet<>());
                                            }

                                            for (int di : d) {
                                                cachedDists.get(c).get(newkey).add(di + distFrom5pEnd.get(key.getFirst()) + distFrom5pEnd.get(key.getSecond()) - distFrom5pEnd.get(newkey.getFirst()) - distFrom5pEnd.get(newkey.getSecond()));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                if (cachedDists.get(c).containsKey(key)) {
                    Set<Integer> dists = cachedDists.get(c).get(key);

                    for (Integer dist : dists) {
                        allDists.get(c).add(dist + distFrom5pEnd.get(key.getFirst()) + distFrom5pEnd.get(key.getSecond()) + fq2.length() - 1);
                    }
                }
            }
        }

        return new PairedEndAlignmentInfo(sai1, sai2, allDists);
    }

    private Set<CortexEdge> filterEdges(Set<CortexEdge> edges, int color) {
        Set<CortexEdge> filteredEdges = new HashSet<>();
        for (CortexEdge e : edges) {
            if (e.getColor() == color) {
                filteredEdges.add(e);
            }
        }

        return filteredEdges;
    }

    private List<String> reverseComplement(List<String> kmers) {
        List<String> revComps = new ArrayList<>();

        for (int i = kmers.size() - 1; i >= 0; i--) {
            String kmer = kmers.get(i);
            revComps.add(SequenceUtils.reverseComplement(kmer));
        }

        return revComps;
    }
}
