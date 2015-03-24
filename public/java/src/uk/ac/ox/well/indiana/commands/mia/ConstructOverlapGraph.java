package uk.ac.ox.well.indiana.commands.mia;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

public class ConstructOverlapGraph extends Module {
    @Argument(fullName="lasTable", shortName="l", doc="Alignment table file (.las.txt)")
    public File LAS_TABLE;

    @Argument(fullName="readMap", shortName="m", doc="Read map")
    public File READ_MAP;

    @Argument(fullName="reads", shortName="r", doc="Reads")
    public FastaSequenceFile READS;

    @Argument(fullName="lengthThreshold", shortName="t", doc="Length threshold")
    public Integer LENGTH_THRESHOLD = 8000;

    @Output
    public File out;

    @Output(fullName="rout", shortName="ro", doc="Reduced graph")
    public File rout;

    private Map<Integer, String> reads;
    private Map<Integer, Interval> readMap;

    private class AlignmentEntry {
        private int readA;
        private int readB;
        private boolean readBIsReversed;
        private int startA;
        private int endA;
        private int startB;
        private int endB;
        private int numDiffs;
        private int numTracePoints;

        public AlignmentEntry(String line) {
            line = line.replaceAll("\\.\\.", " ");
            line = line.replaceAll("^\\s+", "");
            line = line.replaceAll("[\\[\\],\\(\\)x:<]", "");
            line = line.replaceAll("diffs", "");
            line = line.replaceAll("trace pts", "");

            String[] fields = line.split("\\s+");

            readA = Integer.valueOf(fields[0]);
            readB = Integer.valueOf(fields[1]);
            readBIsReversed = fields[2].equals("c");
            startA = Integer.valueOf(fields[3]);
            endA = Integer.valueOf(fields[4]);
            startB = Integer.valueOf(fields[5]);
            endB = Integer.valueOf(fields[6]);
            numDiffs = Integer.valueOf(fields[7]);
            numTracePoints = Integer.valueOf(fields[8]);
        }

        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb
                    .append(readA).append("\t")
                    .append(readB).append("\t")
                    .append(readBIsReversed ? "c" : "n").append("\t")
                    .append("[").append(startA).append("..").append(endA).append("]")
                    .append(" x ")
                    .append("[").append(startB).append("..").append(endB).append("]")
                    .append(" : < ")
                    .append(numDiffs).append(" diffs")
                    .append(" (").append(numTracePoints).append(" trace pts)");

            return sb.toString();
        }
    }

    private class Read {
        private int readId;
        private boolean isReversed;
        private String sequence;

        public Read(int readId, boolean isReversed, String sequence) {
            this.readId = readId;
            this.isReversed = isReversed;
            this.sequence = sequence;
        }

        public int getReadId() { return readId; }

        public String getLabel() { return readId + "_" + (isReversed ? "rc" : "fw"); }

        public int getReadLength() {
            return sequence.length();
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Read read = (Read) o;

            if (isReversed != read.isReversed) return false;
            if (readId != read.readId) return false;
            if (sequence != null ? !sequence.equals(read.sequence) : read.sequence != null) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = readId;
            result = 31 * result + (isReversed ? 1 : 0);
            result = 31 * result + (sequence != null ? sequence.hashCode() : 0);
            return result;
        }
    }

    private enum OverlapType { DOVETAIL, CONTAINMENT, TRANSITIVE, CUT, UNKNOWN }

    private class Overlap extends DefaultEdge {
        private OverlapType ot;

        public Overlap() {
            mark(OverlapType.UNKNOWN);
        }

        public Overlap(OverlapType ot) {
            mark(ot);
        }

        public void mark(OverlapType ot) {
            this.ot = ot;
        }

        public boolean isDovetail() { return this.ot == OverlapType.DOVETAIL; }
        public boolean isContainment() { return this.ot == OverlapType.CONTAINMENT; }
        public boolean isTransitive() { return this.ot == OverlapType.TRANSITIVE; }
    }

    private void loadReads() {
        reads = new HashMap<Integer, String>();

        ReferenceSequence rseq;
        while ((rseq = READS.nextSequence()) != null) {
            String[] names = rseq.getName().split("[\\/\\s]+");
            String seq = new String(rseq.getBases());

            reads.put(Integer.valueOf(names[1]), seq);
        }
    }

    private void loadReadMap() {
        readMap = new HashMap<Integer, Interval>();

        LineReader lr = new LineReader(READ_MAP);
        String line;

        int readNum = 1;
        while ((line = lr.getNextRecord()) != null) {
            line = line.replaceAll("^\\s+", "");
            String[] fields = line.split("\\s+");

            int pos0 = Integer.valueOf(fields[0]);
            int pos1 = Integer.valueOf(fields[1]);

            int start = (pos0 < pos1) ? pos0 : pos1;
            int end = (pos0 < pos1) ? pos1 : pos0;

            Interval interval = new Interval("ref", start, end, pos1 < pos0, String.valueOf(readNum));

            readMap.put(readNum, interval);

            readNum++;
        }
    }

    private void printGraph(DirectedGraph<Read, Overlap> g, PrintStream fout, boolean simplifyGraph) {
        Set<Read> invisibleReads = new HashSet<Read>();

        if (simplifyGraph) {
            for (Overlap e : g.edgeSet()) {
                if (e.isContainment()) {
                    Read vt = g.getEdgeTarget(e);

                    invisibleReads.add(vt);
                }
            }
        }

        fout.println("digraph G {");

        for (Read v : g.vertexSet()) {
            String extra = "";
            if (v.getReadId() == 272) {
                extra = " [style=filled, fillcolor=red]";
            } else if (invisibleReads.contains(v)) {
                extra = " [style=invis]";
            }

            fout.println("    \"" + v.getLabel() + "\"" + extra + ";");
        }

        for (Overlap e : g.edgeSet()) {
            Read vs = g.getEdgeSource(e);
            Read vt = g.getEdgeTarget(e);

            String extra = "";
            if (simplifyGraph && !e.isDovetail()) {
                extra = " [style=invis]";
            }

            fout.println("    \"" + vs.getLabel() + "\" -> \"" + vt.getLabel() + "\"" + extra + ";");
        }

        fout.println("}");
    }

    private Set<Read> getSourceVertices(DirectedGraph<Read, Overlap> g, int level, Read v) {
        Set<Read> sourceVertices = new HashSet<Read>();

        Set<Overlap> inEdges = g.incomingEdgesOf(v);
        for (Overlap e : inEdges) {
            Read sv = g.getEdgeSource(e);
            sourceVertices.add(sv);

            if (level > 0) {
                sourceVertices.addAll(getSourceVertices(g, level - 1, sv));
            }
        }

        return sourceVertices;
    }

    private Set<Read> getTargetVertices(DirectedGraph<Read, Overlap> g, int level, Read v) {
        Set<Read> targetVertices = new HashSet<Read>();

        Set<Overlap> outEdges = g.outgoingEdgesOf(v);
        for (Overlap e : outEdges) {
            Read tv = g.getEdgeTarget(e);
            targetVertices.add(tv);

            if (level > 0) {
                targetVertices.addAll(getTargetVertices(g, level - 1, tv));
            }
        }

        return targetVertices;
    }

    private void markTransitivelyInferribleEdges(DirectedGraph<Read, Overlap> g, int level) {
        for (Read v : g.vertexSet()) {
            Set<Read> sourceVertices = getSourceVertices(g, level, v);
            Set<Read> targetVertices = getTargetVertices(g, level, v);

            for (Read sv : sourceVertices) {
                for (Read tv : targetVertices) {
                    if (g.containsEdge(sv, tv)) {
                        g.getEdge(sv, tv).mark(OverlapType.TRANSITIVE);
                    }
                }
            }
        }
    }

    private void markAmbiguousEdgesToCut(Set<Overlap> edges) {
        int numDovetailEdges = 0;
        for (Overlap e : edges) {
            if (e.isDovetail()) {
                numDovetailEdges++;
            }
        }

        if (numDovetailEdges > 1) {
            for (Overlap e : edges) {
                if (e.isDovetail()) {
                    e.mark(OverlapType.CUT);
                }
            }
        }
    }

    private void markAmbiguousEdgesToCut(DirectedGraph<Read, Overlap> g) {
        for (Read v : g.vertexSet()) {
            markAmbiguousEdgesToCut(g.incomingEdgesOf(v));
            markAmbiguousEdgesToCut(g.outgoingEdgesOf(v));
        }
    }

    private void removeNonDovetailEdges(DirectedGraph<Read, Overlap> g) {
        Set<Overlap> markedEdges = new HashSet<Overlap>();

        for (Overlap e : g.edgeSet()) {
            if (!e.isDovetail()) {
                markedEdges.add(e);
            }
        }

        g.removeAllEdges(markedEdges);
    }

    @Override
    public void execute() {
        loadReads();
        loadReadMap();

        LineReader lr = new LineReader(LAS_TABLE);
        String line;

        List<AlignmentEntry> alignmentEntries = new ArrayList<AlignmentEntry>();
        Map<Pair<Integer, Integer>, Integer> alignmentCounts = new HashMap<Pair<Integer, Integer>, Integer>();
        Set<Integer> containedReads = new HashSet<Integer>();
        while ((line = lr.getNextRecord()) != null) {
            if (line.contains("diffs") && line.contains("trace")) {
                AlignmentEntry ae = new AlignmentEntry(line);
                alignmentEntries.add(ae);

                Pair<Integer, Integer> ap = new Pair<Integer, Integer>(ae.readA, ae.readB);
                if (!alignmentCounts.containsKey(ap)) {
                    alignmentCounts.put(ap, 1);
                } else {
                    alignmentCounts.put(ap, alignmentCounts.get(ap) + 1);
                }

                Read readA = new Read(ae.readA, false, reads.get(ae.readA));
                Read readB = new Read(ae.readB, ae.readBIsReversed, ae.readBIsReversed ? SequenceUtils.reverseComplement(reads.get(ae.readB)) : reads.get(ae.readB));

                if (ae.startA == 0 && ae.endA == readA.getReadLength()) {
                    containedReads.add(ae.readA);
                } else  if (ae.startB == 0 && ae.endB == readB.getReadLength()) {
                    containedReads.add(ae.readB);
                }
            }
        }

        DirectedGraph<Read, Overlap> g = new DefaultDirectedGraph<Read, Overlap>(Overlap.class);

        log.info("Building overlap graph...");
        for (AlignmentEntry ae : alignmentEntries) {
            if ((ae.readA == 114 && ae.readB == 250) || (ae.readB == 114 && ae.readA == 250)) {
                log.info("Hi!");
            }

            Pair<Integer, Integer> ap = new Pair<Integer, Integer>(ae.readA, ae.readB);

            if (alignmentCounts.get(ap) == 1) {
                Read readA = new Read(ae.readA, false, reads.get(ae.readA));
                Read readB = new Read(ae.readB, ae.readBIsReversed, ae.readBIsReversed ? SequenceUtils.reverseComplement(reads.get(ae.readB)) : reads.get(ae.readB));

                Read readARev = new Read(ae.readA, true, SequenceUtils.reverseComplement(readA.sequence));
                Read readBRev = new Read(ae.readB, !ae.readBIsReversed, SequenceUtils.reverseComplement(readB.sequence));

                if (readA.getReadLength() >= LENGTH_THRESHOLD && readB.getReadLength() >= LENGTH_THRESHOLD && !containedReads.contains(ae.readA) && !containedReads.contains(ae.readB)) {
                    //log.info("{} {} {}", ae, readA.getReadLength(), readB.getReadLength());

                    if (ae.startA == 0 && ae.endB == readB.getReadLength()) {
                        // overlap (A extends B)
                        //       ----------->
                        // ----------->

                        g.addVertex(readA);
                        g.addVertex(readB);
                        g.addVertex(readARev);
                        g.addVertex(readBRev);
                        g.addEdge(readB, readA, new Overlap(OverlapType.DOVETAIL));
                        g.addEdge(readARev, readBRev, new Overlap(OverlapType.DOVETAIL));
                    } else if (ae.endA == readA.getReadLength() && ae.startB == 0) {
                        // overlap (B extends A)
                        // ---------->
                        //        ---------->

                        g.addVertex(readA);
                        g.addVertex(readB);
                        g.addVertex(readARev);
                        g.addVertex(readBRev);
                        g.addEdge(readA, readB, new Overlap(OverlapType.DOVETAIL));
                        g.addEdge(readBRev, readARev, new Overlap(OverlapType.DOVETAIL));
                    } else if (ae.startA == 0 && ae.endA == readA.getReadLength()) {
                        // containment (B contains A)
                        //    ---------->
                        // ----------------->

                        g.addVertex(readA);
                        g.addVertex(readB);
                        g.addVertex(readARev);
                        g.addVertex(readBRev);
                        g.addEdge(readB, readA, new Overlap(OverlapType.CONTAINMENT));
                        g.addEdge(readARev, readBRev, new Overlap(OverlapType.CONTAINMENT));
                    } else if (ae.startB == 0 && ae.endB == readB.getReadLength()) {
                        // containment (A contains B)
                        // ----------------->
                        //    ---------->

                        g.addVertex(readA);
                        g.addVertex(readB);
                        g.addVertex(readARev);
                        g.addVertex(readBRev);
                        g.addEdge(readA, readB, new Overlap(OverlapType.CONTAINMENT));
                        g.addEdge(readBRev, readARev, new Overlap(OverlapType.CONTAINMENT));
                    } else {
                        // chimeric?
                        // xx-------------xx>
                        // xx-------------xx>
                    }
                }
            }
        }

        log.info("  stats (before layout)");
        log.info("  - nodes: {}", g.vertexSet().size());
        log.info("  - edges: {}", g.edgeSet().size());

        try {
            printGraph(g, new PrintStream(out), false);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        markTransitivelyInferribleEdges(g, 3);
        markAmbiguousEdgesToCut(g);
        //removeNonDovetailEdges(g);

        log.info("  stats (after layout)");
        log.info("  - nodes: {}", g.vertexSet().size());
        log.info("  - edges: {}", g.edgeSet().size());

        try {
            printGraph(g, new PrintStream(rout), true);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}
