package uk.ac.ox.well.indiana.commands.caller.call;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.traverse.TopologicalOrderIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;

import java.util.*;

/**
 * Created by kiran on 23/06/2017.
 */
public class MergeContigs extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="sequences", shortName="s", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="drafts", shortName="d", doc="Drafts")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Argument(fullName="verify", shortName="v", doc="Verify")
    public FastaSequenceFile VERIFY;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(PARENTS);
        List<Integer> recruitColors = GRAPH.getColorsForSampleNames(new ArrayList<>(LOOKUPS.keySet()));

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(childColor)
                .joiningColors(parentColors)
                .recruitmentColors(recruitColors)
                .graph(GRAPH)
                .make();

        ReferenceSequence vseq;
        Set<CortexKmer> validatedKmers = new HashSet<>();
        while ((vseq = VERIFY.nextSequence()) != null) {
            for (int i = 0; i <= vseq.length() - GRAPH.getKmerSize(); i++) {
                CortexKmer ck = new CortexKmer(vseq.getBaseString().substring(i, i + GRAPH.getKmerSize()));
                validatedKmers.add(ck);
            }
        }

        DirectedGraph<ReferenceSequence, LabeledEdge> g = loadContigs();
        Map<ReferenceSequence, Integer> mergeable = new HashMap<>();

        for (ReferenceSequence rseq : g.vertexSet()) {
            if (!rseq.getName().contains("boundary")) {
                for (int i = 0; i <= rseq.length() - GRAPH.getKmerSize(); i++) {
                    CortexKmer ck = new CortexKmer(rseq.getBaseString().substring(i, i + GRAPH.getKmerSize()));
                    if (validatedKmers.contains(ck)) {
                        //mergeable.add(rseq);
                        ContainerUtils.increment(mergeable, rseq);
                    }
                }
            }
        }

        log.info("Mergeable:");
        for (ReferenceSequence rseq : mergeable.keySet()) {
            if (mergeable.get(rseq) > 20) {
                log.info("  {} {}", mergeable.get(rseq), rseq);
            }
        }
        log.info("");

        log.info("Processing contigs:");

        Set<ReferenceSequence> toRemove = new HashSet<>();
        for (ReferenceSequence rseq : g.vertexSet()) {
            if (rseq.getContigIndex() > -1) {
                if (rseq.getName().contains("contig32") || rseq.getName().contains("contig97")) {
                    log.info("  {}", rseq);

                    for (int i = 0; i <= rseq.length() - GRAPH.getKmerSize(); i++) {
                        String sk = rseq.getBaseString().substring(i, i + GRAPH.getKmerSize());
                        CortexKmer ck = new CortexKmer(sk);

                        if (validatedKmers.contains(ck)) {
                            log.info("    - {}/{} {} {}", i, rseq.length() - GRAPH.getKmerSize(), sk, validatedKmers.contains(ck));
                        }
                    }
                }

                String adjRev = extend(e, g, rseq, false);
                String adjFwd = extend(e, g, rseq, true);

                log.info("  {}", rseq.getName());
                log.info("    - {}", adjRev);
                log.info("    - {}", adjFwd);
            } else {
                toRemove.add(rseq);
            }
        }

        g.removeAllVertices(toRemove);
    }

    private String extend(TraversalEngine e, DirectedGraph<ReferenceSequence, LabeledEdge> g, ReferenceSequence rseq, boolean goForward) {
        Map<String, Interval> mostRecentConfidentInterval = new HashMap<>();
        Set<String> acceptableContigs = new HashSet<>();

        int start = !goForward ? 0 : rseq.length() - GRAPH.getKmerSize();
        int end   = !goForward ? rseq.length() - GRAPH.getKmerSize() : 0;
        int inc   = !goForward ? 1 : -1;

        for (int i = start; !goForward ? i <= end : i >= end; i += inc) {
            String sk = rseq.getBaseString().substring(i, i + GRAPH.getKmerSize());

            for (String background : LOOKUPS.keySet()) {
                Set<Interval> intervals = LOOKUPS.get(background).findKmer(sk);

                if (intervals.size() == 1 && !mostRecentConfidentInterval.containsKey(background)) {
                    Interval it = intervals.iterator().next();
                    mostRecentConfidentInterval.put(background, it);
                    acceptableContigs.add(it.getContig());
                }

                if (mostRecentConfidentInterval.size() == LOOKUPS.size()) {
                    break;
                }
            }
        }

        if (mostRecentConfidentInterval.size() > 0) {
            List<ReferenceSequence> arseqs = !goForward ? Graphs.predecessorListOf(g, rseq) : Graphs.successorListOf(g, rseq);

            if (arseqs.size() == 1) {
                ReferenceSequence arseq = arseqs.get(0);

                if ((!goForward && Graphs.vertexHasPredecessors(g, arseq)) || (goForward && Graphs.vertexHasSuccessors(g, arseq))) {
                    return null;
                }

                String sk = arseq.getBaseString();
                List<String> gap = new ArrayList<>();

                do {
                    Set<CortexVertex> avs = !goForward ? e.getPrevVertices(sk) : e.getNextVertices(sk);

                    sk = null;

                    if (avs.size() == 1) {
                        CortexVertex av = avs.iterator().next();
                        sk = av.getSk();
                    } else if (avs.size() > 1) {
                        Set<String> adj = new HashSet<>();
                        for (CortexVertex av : avs) {
                            for (String background : LOOKUPS.keySet()) {
                                Set<Interval> its = LOOKUPS.get(background).findKmer(av.getSk());

                                if (its.size() == 1) {
                                    Interval it = its.iterator().next();
                                    if (acceptableContigs.contains(it.getContig())) {
                                        adj.add(av.getSk());
                                    }
                                }
                            }
                        }

                        if (adj.size() == 1) {
                            sk = adj.iterator().next();
                        }
                    }

                    if (sk != null) {
                        ReferenceSequence boundary = new ReferenceSequence("boundary", -1, sk.getBytes());
                        if (g.containsVertex(boundary)) {
                            StringBuilder sb = new StringBuilder();
                            for (String s : gap) {
                                if (sb.length() == 0) {
                                    sb.append(s.toLowerCase());
                                } else {
                                    sb.append(s.substring(s.length() - 1, s.length()).toLowerCase());
                                }
                            }

                            List<ReferenceSequence> arseqs2 = !goForward ? Graphs.predecessorListOf(g, rseq) : Graphs.successorListOf(g, rseq);
                            ReferenceSequence aseq = arseqs2.get(0);

                            if (!goForward) {
                                g.addEdge(aseq, rseq, new LabeledEdge(sb.toString()));

                                log.info("  joined");
                                log.info("    {}", aseq.getName());
                                log.info("    {}", rseq.getName());
                            } else {
                                g.addEdge(rseq, aseq, new LabeledEdge(sb.toString()));

                                log.info("  joined");
                                log.info("    {}", rseq.getName());
                                log.info("    {}", aseq.getName());
                            }

                            return aseq.getName();
                        } else {
                            if (!goForward) {
                                gap.add(0, sk);
                            } else {
                                gap.add(sk);
                            }
                        }
                    }
                } while (sk != null);
            }
        }

        return null;
    }

    private class LabeledEdge extends DefaultEdge {
        private String label;

        public LabeledEdge() {}
        public LabeledEdge(String label) { this.label = label; }

        public void setLabel(String label) { this.label = label; }
        public String getLabel() { return label; }
    }

    private DirectedGraph<ReferenceSequence, LabeledEdge> loadContigs() {
        DirectedGraph<ReferenceSequence, LabeledEdge> g = new DefaultDirectedGraph<>(LabeledEdge.class);

        ReferenceSequence fwseq;
        while ((fwseq = CONTIGS.nextSequence()) != null) {
            String fw = fwseq.getBaseString();
            String rc = SequenceUtils.reverseComplement(fw);
            ReferenceSequence rcseq   = new ReferenceSequence(fwseq.getName(), fwseq.getContigIndex(), rc.getBytes());

            String fwFirstKmer = fwseq.getBaseString().substring(0, GRAPH.getKmerSize());
            ReferenceSequence fwFirst = new ReferenceSequence(fwFirstKmer, -1, fwFirstKmer.getBytes());
            String fwLastKmer = fwseq.getBaseString().substring(fwseq.length() - GRAPH.getKmerSize(), fwseq.length());
            ReferenceSequence fwLast  = new ReferenceSequence(fwLastKmer, -1, fwLastKmer.getBytes());

            String rcFirstKmer = rcseq.getBaseString().substring(0, GRAPH.getKmerSize());
            ReferenceSequence rcFirst = new ReferenceSequence(rcFirstKmer, -1, rcFirstKmer.getBytes());
            String rcLastKmer = rcseq.getBaseString().substring(rcseq.length() - GRAPH.getKmerSize(), rcseq.length());
            ReferenceSequence rcLast  = new ReferenceSequence(rcLastKmer, -1, rcLastKmer.getBytes());

            g.addVertex(fwFirst);
            g.addVertex(fwseq);
            g.addVertex(fwLast);
            g.addEdge(fwFirst, fwseq, new LabeledEdge());
            g.addEdge(fwseq, fwLast, new LabeledEdge());

            g.addVertex(rcFirst);
            g.addVertex(rcseq);
            g.addVertex(rcLast);
            g.addEdge(rcFirst, rcseq, new LabeledEdge());
            g.addEdge(rcseq, rcLast, new LabeledEdge());
        }

        return g;
    }
}
