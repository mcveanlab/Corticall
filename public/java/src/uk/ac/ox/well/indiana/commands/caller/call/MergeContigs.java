package uk.ac.ox.well.indiana.commands.caller.call;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;

import java.io.PrintStream;
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

    @Output
    public PrintStream out;

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

        DirectedGraph<Contig, LabeledEdge> g = loadContigs();

        mergeContigs(e, g);

        emitContigs(g);
    }

    private void mergeContigs(TraversalEngine e, DirectedGraph<Contig, LabeledEdge> g) {
        log.info("Processing contigs:");

        Set<Contig> toRemove = new HashSet<>();
        for (Contig rseq : g.vertexSet()) {
            if (rseq.getIndex() > -1) {
                log.info("  {}", rseq.getName().split("\\s+")[0]);

                String adjRev = extend(e, g, rseq, false);
                if (adjRev != null) {
                    log.info("    - joined prev {}", adjRev);
                }

                String adjFwd = extend(e, g, rseq, true);
                if (adjFwd != null) {
                    log.info("    - joined next {}", adjFwd);
                }
            } else {
                toRemove.add(rseq);
            }
        }
        g.removeAllVertices(toRemove);
    }

    private void emitContigs(DirectedGraph<Contig, LabeledEdge> g) {
        Set<String> seen = new HashSet<>();
        for (Contig rseq : g.vertexSet()) {
            if (!seen.contains(rseq.getName().split("\\s+")[0])) {
                Set<String> names = new TreeSet<>();
                List<String> scaffold = new ArrayList<>();
                scaffold.add(rseq.getSequence());

                names.add(rseq.getName().split("\\s+")[0]);

                Contig cur = rseq;
                while (Graphs.vertexHasPredecessors(g, cur)) {
                    Contig pre = Graphs.predecessorListOf(g, cur).get(0);
                    //Contig gap = new Contig(g.getEdge(pre, cur).getLabel());

                    names.add(pre.getName().split("\\s+")[0]);

                    int overlap = g.getEdge(pre, cur).getOverlap();
                    if (overlap > 0) {
                        //pre = new Contig(pre.getName(), pre.getSequence().substring(0, pre.length() - overlap), pre.getIndex());

                        scaffold.add(0, pre.getSequence().substring(0, pre.length() - overlap));
                    } else {
                        scaffold.add(0, pre.getSequence());
                    }

                    //scaffold.add(0, gap);
                    //scaffold.add(0, pre);

                    cur = pre;
                }

                cur = rseq;
                while (Graphs.vertexHasSuccessors(g, cur)) {
                    Contig nxt = Graphs.successorListOf(g, cur).get(0);
                    //Contig gap = new Contig(g.getEdge(cur, nxt).getLabel());

                    names.add(nxt.getName().split("\\s+")[0]);

                    int overlap = g.getEdge(cur, nxt).getOverlap();
                    if (overlap > 0) {
                        //nxt = new Contig(nxt.getName(), nxt.getSequence().substring(overlap, nxt.length()), nxt.getIndex());

                        scaffold.add(nxt.getSequence().substring(overlap, nxt.length()));
                    } else {
                        scaffold.add(nxt.getSequence());
                    }

                    //scaffold.add(gap);
                    //scaffold.add(nxt);

                    cur = nxt;
                }

                StringBuilder sb = new StringBuilder();
                String name = Joiner.on("_").join(names);
                for (String c : scaffold) {
                    //String[] pieces = c.getName().split("\\s+");
                    //name.append(pieces[0]);

                    sb.append(c);
                }

                out.println(">" + name);
                out.println(sb);

                for (String n : names) {
                    seen.add(n);
                }
            }
        }
    }

    private String extend(TraversalEngine e, DirectedGraph<Contig, LabeledEdge> g, Contig rseq, boolean goForward) {
        Map<String, Interval> mostRecentConfidentInterval = new HashMap<>();
        Set<String> acceptableContigs = new HashSet<>();

        int start = !goForward ? 0 : rseq.length() - GRAPH.getKmerSize();
        int end   = !goForward ? rseq.length() - GRAPH.getKmerSize() : 0;
        int inc   = !goForward ? 1 : -1;

        for (int i = start; !goForward ? i <= end : i >= end; i += inc) {
            String sk = rseq.getSequence().substring(i, i + GRAPH.getKmerSize());

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
            List<Contig> arseqs = !goForward ? Graphs.predecessorListOf(g, rseq) : Graphs.successorListOf(g, rseq);

            if (arseqs.size() > 1) {
                for (Contig arseq : arseqs) {
                    if (arseq.getIndex() > -1 && !goForward && Graphs.vertexHasPredecessors(g, arseq)) {
                        return arseq.getName();
                    }

                    if (arseq.getIndex() > -1 && goForward && Graphs.vertexHasSuccessors(g, arseq)) {
                        return arseq.getName();
                    }
                }
            } else if (arseqs.size() == 1) {
                Contig arseq = arseqs.get(0);

                if ((!goForward && Graphs.vertexHasPredecessors(g, arseq)) || (goForward && Graphs.vertexHasSuccessors(g, arseq))) {
                    return null;
                }

                String sk = arseq.getSequence();
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
                        Contig boundary = new Contig(sk);
                        if (g.containsVertex(boundary)) {
                            StringBuilder sb = new StringBuilder();
                            for (String s : gap) {
                                if (sb.length() == 0) {
                                    sb.append(s.toLowerCase());
                                } else {
                                    sb.append(s.substring(s.length() - 1, s.length()).toLowerCase());
                                }
                            }

                            int intersectionLength = -1;
                            for (String background : LOOKUPS.keySet()) {
                                Set<Interval> its = LOOKUPS.get(background).findKmer(sk);
                                if (its.size() == 1 && mostRecentConfidentInterval.containsKey(background)) {
                                    Interval it1 = mostRecentConfidentInterval.get(background);
                                    Interval it2 = its.iterator().next();
                                    int newIntersectionLength = it1.getIntersectionLength(it2);

                                    if (intersectionLength == -1 || newIntersectionLength > 0) {
                                        intersectionLength = newIntersectionLength;
                                    }
                                }
                            }

                            // Testing
                            if (rseq.getName().contains("contig58")) {
                                Contig aseq = Graphs.successorListOf(g, boundary).get(0);

                                for (int i = 0; i < Math.max(rseq.length(), aseq.length()) - GRAPH.getKmerSize(); i++) {
                                    String skr = (i + GRAPH.getKmerSize() <= rseq.length()) ? rseq.getSequence().substring(i, i + GRAPH.getKmerSize()) : null;
                                    String ska = (i + GRAPH.getKmerSize() <= aseq.length()) ? aseq.getSequence().substring(i, i + GRAPH.getKmerSize()) : null;

                                    Set<Interval> isr = new TreeSet<>();
                                    if (skr != null) {
                                        for (String background : LOOKUPS.keySet()) {
                                            isr.addAll(LOOKUPS.get(background).findKmer(skr));
                                        }
                                    }

                                    Set<Interval> isa = new TreeSet<>();
                                    if (ska != null) {
                                        for (String background : LOOKUPS.keySet()) {
                                            isa.addAll(LOOKUPS.get(background).findKmer(ska));
                                        }
                                    }

                                    log.info("{} {} {} {} {}", i, skr, ska, isr, isa);
                                }
                            }
                            // Testing

                            log.info("Intersection: {}", intersectionLength);

                            List<Contig> arseqs2 = !goForward ? Graphs.predecessorListOf(g, boundary) : Graphs.successorListOf(g, boundary);
                            Contig aseq = arseqs2.get(0);
                            //String gapString = sb.toString();

                            if (!goForward) {
                                //g.addEdge(aseq, rseq, new LabeledEdge(gapString.substring(GRAPH.getKmerSize() - 1, gapString.length())));
                                g.addEdge(aseq, rseq, new LabeledEdge(intersectionLength));

                                log.info("  joined");
                                log.info("    {}", aseq.getName());
                                log.info("    {}", rseq.getName());
                            } else {
                                //g.addEdge(rseq, aseq, new LabeledEdge(gapString.substring(0, gapString.length() - GRAPH.getKmerSize() + 1)));
                                g.addEdge(rseq, aseq, new LabeledEdge(intersectionLength));

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

    private class Contig {
        private String name;
        private String sequence;
        private int index;

        public Contig(ReferenceSequence rseq) {
            name = rseq.getName();
            sequence = rseq.getBaseString();
            index = rseq.getContigIndex();
        }

        public Contig(String seq) {
            name = seq;
            sequence = seq;
            index = -1;
        }

        public Contig(String name, String seq, int index) {
            this.name = name;
            this.sequence = seq;
            this.index = index;
        }

        public String getName() { return name; }
        public void setName(String name) { this.name = name; }

        public String getSequence() { return sequence; }
        public void setSequence(String sequence) { this.sequence = sequence; }

        public int getIndex() { return index; }
        public void setIndex(int index) { this.index = index; }

        public int length() { return sequence.length(); }

        @Override
        public String toString() {
            return "Contig{" +
                    "name='" + name + '\'' +
                    ", sequence='" + (sequence.length() > 100 ? sequence.substring(0, 100) : sequence) + '\'' +
                    ", index=" + index +
                    '}';
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Contig contig = (Contig) o;

            if (index != contig.index) return false;
            if (name != null ? !name.equals(contig.name) : contig.name != null) return false;
            return sequence != null ? sequence.equals(contig.sequence) : contig.sequence == null;

        }

        @Override
        public int hashCode() {
            int result = name != null ? name.hashCode() : 0;
            result = 31 * result + (sequence != null ? sequence.hashCode() : 0);
            result = 31 * result + index;
            return result;
        }
    }

    private class LabeledEdge extends DefaultEdge {
        private String label = "";
        private int overlap = 0;

        public LabeledEdge() {}
        public LabeledEdge(String label) { this.label = label; }
        public LabeledEdge(int overlap) { this.overlap = overlap; }
        public LabeledEdge(String label, int overlap) { this.label = label; this.overlap = overlap; }

        public void setLabel(String label) { this.label = label; }
        public String getLabel() { return label; }

        public void setOverlap(int overlap) { this.overlap = overlap; }
        public int getOverlap() { return overlap; }
    }

    private DirectedGraph<Contig, LabeledEdge> loadContigs() {
        DirectedGraph<Contig, LabeledEdge> g = new DefaultDirectedGraph<>(LabeledEdge.class);

        ReferenceSequence seq;
        while ((seq = CONTIGS.nextSequence()) != null) {
            Contig fwseq = new Contig(seq);
            Contig rcseq = new Contig(new ReferenceSequence(seq.getName(), seq.getContigIndex(), SequenceUtils.reverseComplement(seq.getBases())));

            Contig fwFirst = new Contig(fwseq.getSequence().substring(0, GRAPH.getKmerSize()));
            Contig fwLast  = new Contig(fwseq.getSequence().substring(fwseq.length() - GRAPH.getKmerSize(), fwseq.length()));

            Contig rcFirst = new Contig(rcseq.getSequence().substring(0, GRAPH.getKmerSize()));
            Contig rcLast  = new Contig(rcseq.getSequence().substring(rcseq.length() - GRAPH.getKmerSize(), rcseq.length()));

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
