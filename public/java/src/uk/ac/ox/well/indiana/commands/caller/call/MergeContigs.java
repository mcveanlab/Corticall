package uk.ac.ox.well.indiana.commands.caller.call;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.traverse.TopologicalOrderIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
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

        //DirectedGraph<ReferenceSequence, String> q = loadContigs();
        DirectedGraph<ReferenceSequence, String> g = loadContigs();

        log.info("Processing contigs:");

        Set<ReferenceSequence> toRemove = new HashSet<>();
        for (ReferenceSequence rseq : g.vertexSet()) {
            if (rseq.getContigIndex() > -1) {
                //log.info("  {}", rseq);

                extend(e, g, rseq, false);
                extend(e, g, rseq, true);
            } else {
                toRemove.add(rseq);
            }
        }

        g.removeAllVertices(toRemove);

        TopologicalOrderIterator<ReferenceSequence, String> toi = new TopologicalOrderIterator<ReferenceSequence, String>(g);
        while (toi.hasNext()) {
            log.info("{}", toi.next());
        }
    }

    private void extend(TraversalEngine e, DirectedGraph<ReferenceSequence, String> g, ReferenceSequence rseq, boolean goForward) {
        Map<String, Interval> mostRecentConfidentInterval = new HashMap<>();
        Set<String> acceptableContigs = new HashSet<>();

        int start = !goForward ? 0 : rseq.length() - GRAPH.getKmerSize();
        int end   = !goForward ? rseq.length() - GRAPH.getKmerSize() : 0;
        int inc   = !goForward ? 1 : -1;

        //log.info("{} {} {} {}", goForward, start, end, inc);

        //for (int i = rseq.length() - GRAPH.getKmerSize(); i >= 0; i--) {
        for (int i = start; !goForward ? i <= end : i >= end; i += inc) {
            //log.info("  {} {}", i, i + GRAPH.getKmerSize());

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

                if ((!goForward && Graphs.vertexHasPredecessors(g, rseq)) || (goForward && Graphs.vertexHasSuccessors(g, rseq))) {
                    return;
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
                        for (CortexVertex av : avs) {
                            for (String background : LOOKUPS.keySet()) {
                                Set<Interval> its = LOOKUPS.get(background).findKmer(av.getSk());

                                if (its.size() == 1) {
                                    Interval it = its.iterator().next();
                                    if (acceptableContigs.contains(it.getContig())) {
                                        sk = av.getSk();
                                    }
                                }
                            }
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

                            //g.addEdge(rseq, boundary, sb.toString());

                            List<ReferenceSequence> arseqs2 = !goForward ? Graphs.predecessorListOf(g, rseq) : Graphs.successorListOf(g, rseq);
                            ReferenceSequence aseq = arseqs2.get(0);

                            if (!goForward) {
                                g.addEdge(aseq, rseq, sb.toString());

                                log.info("  joined");
                                log.info("    {}", aseq.getName());
                                log.info("    {}", rseq.getName());
                            } else {
                                g.addEdge(rseq, aseq, sb.toString());

                                log.info("  joined");
                                log.info("    {}", rseq.getName());
                                log.info("    {}", aseq.getName());
                            }

                            sk = null;
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
    }

    private DirectedGraph<ReferenceSequence, String> loadContigs() {
        DirectedGraph<ReferenceSequence, String> g = new DefaultDirectedGraph<>(String.class);

        ReferenceSequence fwseq;
        while ((fwseq = CONTIGS.nextSequence()) != null) {
            String fw = fwseq.getBaseString();
            String rc = SequenceUtils.reverseComplement(fw);

            ReferenceSequence fwFirst = new ReferenceSequence("boundary", -1, fwseq.getBaseString().substring(0, GRAPH.getKmerSize()).getBytes());
            ReferenceSequence fwLast  = new ReferenceSequence("boundary", -1, fwseq.getBaseString().substring(fwseq.length() - GRAPH.getKmerSize(), fwseq.length()).getBytes());

            ReferenceSequence rcseq   = new ReferenceSequence(fwseq.getName(), fwseq.getContigIndex(), rc.getBytes());
            ReferenceSequence rcFirst = new ReferenceSequence("boundary", -1, rcseq.getBaseString().substring(0, GRAPH.getKmerSize()).getBytes());
            ReferenceSequence rcLast  = new ReferenceSequence("boundary", -1, rcseq.getBaseString().substring(rcseq.length() - GRAPH.getKmerSize(), rcseq.length()).getBytes());

            g.addVertex(fwFirst);
            g.addVertex(fwseq);
            g.addVertex(fwLast);
            g.addEdge(fwFirst, fwseq);
            g.addEdge(fwseq, fwLast);

            g.addVertex(rcFirst);
            g.addVertex(rcseq);
            g.addVertex(rcLast);
            g.addEdge(rcFirst, rcseq);
            g.addEdge(rcseq, rcLast);
        }

        return g;
    }
}
