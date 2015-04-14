package uk.ac.ox.well.indiana.commands.mia;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
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

public class ConstructChunkGraph extends Module {
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
        private boolean readAIsReversed;
        private int startA;
        private int endA;
        private int lengthA;

        private int readB;
        private boolean readBIsReversed;
        private int startB;
        private int endB;
        private int lengthB;

        private int numDiffs;
        private int numTracePoints;

        public AlignmentEntry() {}

        public AlignmentEntry(String line) {
            line = line.replaceAll("\\.\\.", " ");
            line = line.replaceAll("^\\s+", "");
            line = line.replaceAll("[\\[\\],\\(\\)x:<]", "");
            line = line.replaceAll("diffs", "");
            line = line.replaceAll("trace pts", "");

            String[] fields = line.split("\\s+");

            readA = Integer.valueOf(fields[0]);
            startA = Integer.valueOf(fields[3]);
            endA = Integer.valueOf(fields[4]);
            readAIsReversed = false;
            lengthA = 0;

            readB = Integer.valueOf(fields[1]);
            readBIsReversed = fields[2].equals("c");
            startB = Integer.valueOf(fields[5]);
            endB = Integer.valueOf(fields[6]);
            lengthB = 0;

            numDiffs = Integer.valueOf(fields[7]);
            numTracePoints = Integer.valueOf(fields[8]);
        }

        public void setLengthA(int lengthA) { this.lengthA = lengthA; }
        public void setLengthB(int lengthB) { this.lengthB = lengthB; }

        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb
                    .append(readA).append("\t")
                    .append(readB).append("\t")
                    .append(readAIsReversed ? "c" : "n").append("\t")
                    .append(readBIsReversed ? "c" : "n").append("\t")
                    .append("[").append(startA).append("..").append(endA).append("] (").append(lengthA).append(") ")
                    .append(" x ")
                    .append("[").append(startB).append("..").append(endB).append("] (").append(lengthB).append(") ")
                    .append(" : < ")
                    .append(numDiffs).append(" diffs")
                    .append(" (").append(numTracePoints).append(" trace pts)");

            return sb.toString();
        }

        public AlignmentEntry flip() {
            AlignmentEntry aeFlipped = new AlignmentEntry();

            aeFlipped.readA = readB;
            aeFlipped.readAIsReversed = !readBIsReversed;
            aeFlipped.lengthA = lengthB;
            aeFlipped.startA = lengthB - endB;
            aeFlipped.endA = lengthB - startB;

            aeFlipped.readB = readA;
            aeFlipped.readBIsReversed = !readAIsReversed;
            aeFlipped.lengthB = lengthA;
            aeFlipped.startB = lengthA - endA;
            aeFlipped.endB = lengthA - startA;

            aeFlipped.numDiffs = numDiffs;
            aeFlipped.numTracePoints = numTracePoints;

            return aeFlipped;
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

        public int getReadLength() { return sequence.length(); }

        public String getReadSequence() { return sequence; }

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

        @Override
        public String toString() {
            return "Read{" +
                    "readId=" + readId +
                    ", isReversed=" + isReversed +
                    '}';
        }
    }

    private enum OverlapType { DOVETAIL, CONTAINMENT, TRANSITIVE, CUT, UNKNOWN }

    private class Overlap extends DefaultEdge {
        private OverlapType ot;
        private AlignmentEntry ae;

        public Overlap() { mark(OverlapType.UNKNOWN); }
        public Overlap(OverlapType ot, AlignmentEntry ae) {
            mark(ot);
            this.ae = ae;
        }

        public void mark(OverlapType ot) { this.ot = ot; }

        public boolean isDovetail() { return this.ot == OverlapType.DOVETAIL; }
        public boolean isContainment() { return this.ot == OverlapType.CONTAINMENT; }
        public boolean isTransitive() { return this.ot == OverlapType.TRANSITIVE; }

        public AlignmentEntry getAlignment() {
            return ae;
        }
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
            if (invisibleReads.contains(v)) {
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

    private DirectedGraph<Read, Overlap> constructOverlapGraph() {
        loadReads();
        loadReadMap();

        LineReader lr = new LineReader(LAS_TABLE);
        String line;

        Map<Pair<Integer, Integer>, AlignmentEntry> alignmentEntryMap = new HashMap<Pair<Integer, Integer>, AlignmentEntry>();
        Map<Pair<Integer, Integer>, Integer> alignmentCounts = new HashMap<Pair<Integer, Integer>, Integer>();
        Set<Integer> containedReads = new HashSet<Integer>();

        while ((line = lr.getNextRecord()) != null) {
            if (line.contains("diffs") && line.contains("trace")) {
                AlignmentEntry ae = new AlignmentEntry(line);
                Pair<Integer, Integer> ap = new Pair<Integer, Integer>(ae.readA, ae.readB);

                Read readA = new Read(ae.readA, false, reads.get(ae.readA));
                Read readB = new Read(ae.readB, ae.readBIsReversed, ae.readBIsReversed ? SequenceUtils.reverseComplement(reads.get(ae.readB)) : reads.get(ae.readB));

                ae.setLengthA(readA.getReadLength());
                ae.setLengthB(readB.getReadLength());

                if (!alignmentEntryMap.containsKey(ap)) {
                    alignmentEntryMap.put(ap, ae);
                }

                if (!alignmentCounts.containsKey(ap)) {
                    alignmentCounts.put(ap, 1);
                } else {
                    alignmentCounts.put(ap, alignmentCounts.get(ap) + 1);
                }

                if (ae.startA == 0 && ae.endA == readA.getReadLength()) {
                    containedReads.add(ae.readA);
                } else  if (ae.startB == 0 && ae.endB == readB.getReadLength()) {
                    containedReads.add(ae.readB);
                }
            }
        }

        DirectedGraph<Read, Overlap> g = new DefaultDirectedGraph<Read, Overlap>(Overlap.class);

        for (Pair<Integer, Integer> ap : alignmentEntryMap.keySet()) {
            Pair<Integer, Integer> apc = new Pair<Integer, Integer>(ap.getSecond(), ap.getFirst());

            if (alignmentCounts.get(ap) == 1) {
                AlignmentEntry ae = alignmentEntryMap.get(ap);
                Read readA = new Read(ae.readA, false, reads.get(ae.readA).toUpperCase());
                Read readB = new Read(ae.readB, ae.readBIsReversed, ae.readBIsReversed ? SequenceUtils.reverseComplement(reads.get(ae.readB)) : reads.get(ae.readB).toUpperCase());

                AlignmentEntry aeRev = ae.flip();
                Read readARev = new Read(ae.readA, true, SequenceUtils.reverseComplement(readA.sequence));
                Read readBRev = new Read(ae.readB, !ae.readBIsReversed, SequenceUtils.reverseComplement(readB.sequence));

                log.info("aeFw: {}", ae);
                log.info("aeRc: {}", aeRev);

                log.info("{}", readA.getReadSequence().substring(ae.startA, ae.endA));
                log.info("{}", readB.getReadSequence().substring(ae.startB, ae.endB));

                log.info("{}", readBRev.getReadSequence().substring(aeRev.startA, aeRev.endA));
                log.info("{}", readARev.getReadSequence().substring(aeRev.startB, aeRev.endB));

                if (readA.getReadLength() >= LENGTH_THRESHOLD && readB.getReadLength() >= LENGTH_THRESHOLD && !containedReads.contains(ae.readA) && !containedReads.contains(ae.readB)) {
                    if (ae.startA == 0 && ae.endB == readB.getReadLength()) {
                        // overlap (A extends B)
                        //       ----------->
                        // ----------->

                        g.addVertex(readA);
                        g.addVertex(readB);
                        g.addVertex(readARev);
                        g.addVertex(readBRev);
                        g.addEdge(readB, readA, new Overlap(OverlapType.DOVETAIL, ae));
                        g.addEdge(readARev, readBRev, new Overlap(OverlapType.DOVETAIL, aeRev));
                    } else if (ae.endA == readA.getReadLength() && ae.startB == 0) {
                        // overlap (B extends A)
                        // ---------->
                        //        ---------->

                        g.addVertex(readA);
                        g.addVertex(readB);
                        g.addVertex(readARev);
                        g.addVertex(readBRev);
                        g.addEdge(readA, readB, new Overlap(OverlapType.DOVETAIL, ae));
                        g.addEdge(readBRev, readARev, new Overlap(OverlapType.DOVETAIL, aeRev));
                    } else if (ae.startA == 0 && ae.endA == readA.getReadLength()) {
                        // containment (B contains A)
                        //    ---------->
                        // ----------------->

                        g.addVertex(readA);
                        g.addVertex(readB);
                        g.addVertex(readARev);
                        g.addVertex(readBRev);
                        g.addEdge(readB, readA, new Overlap(OverlapType.CONTAINMENT, ae));
                        g.addEdge(readARev, readBRev, new Overlap(OverlapType.CONTAINMENT, aeRev));
                    } else if (ae.startB == 0 && ae.endB == readB.getReadLength()) {
                        // containment (A contains B)
                        // ----------------->
                        //    ---------->

                        g.addVertex(readA);
                        g.addVertex(readB);
                        g.addVertex(readARev);
                        g.addVertex(readBRev);
                        g.addEdge(readA, readB, new Overlap(OverlapType.CONTAINMENT, ae));
                        g.addEdge(readBRev, readARev, new Overlap(OverlapType.CONTAINMENT, aeRev));
                    } else {
                        // chimeric?
                        // xx-------------xx>
                        // xx-------------xx>
                    }
                }
            }
        }

        return g;
    }

    private int visibleOutDegreeOf(DirectedGraph<Read, Overlap> g, Read v) {
        Set<Overlap> edges = g.outgoingEdgesOf(v);

        int visibleOutDegree = 0;
        for (Overlap e : edges) {
            if (e.isDovetail()) {
                visibleOutDegree++;
            }
        }

        return visibleOutDegree;
    }

    private Set<Overlap> visibleOutgoingEdgesOf(DirectedGraph<Read, Overlap> g, Read v) {
        Set<Overlap> es = new HashSet<Overlap>();

        for (Overlap e : g.outgoingEdgesOf(v)) {
            if (e.isDovetail()) {
                es.add(e);
            }
        }

        return es;
    }

    private Set<Overlap> visibleIncomingEdgesOf(DirectedGraph<Read, Overlap> g, Read v) {
        Set<Overlap> es = new HashSet<Overlap>();

        for (Overlap e : g.incomingEdgesOf(v)) {
            if (e.isDovetail()) {
                es.add(e);
            }
        }

        return es;
    }

    private int visibleInDegreeOf(DirectedGraph<Read, Overlap> g, Read v) {
        Set<Overlap> edges = g.incomingEdgesOf(v);

        int visibleInDegree = 0;
        for (Overlap e : edges) {
            if (e.isDovetail()) {
                visibleInDegree++;
            }
        }

        return visibleInDegree;
    }

    private Set<DirectedGraph<Read, Overlap>> getChunks(DirectedGraph<Read, Overlap> g) {
        Set<Read> startPoints = new HashSet<Read>();

        for (Read v : g.vertexSet()) {
            if (visibleInDegreeOf(g, v) == 0 && visibleOutDegreeOf(g, v) == 1) {
                startPoints.add(v);
            }
        }

        Set<DirectedGraph<Read, Overlap>> chunks = new HashSet<DirectedGraph<Read, Overlap>>();

        for (Read startPoint : startPoints) {
            DirectedGraph<Read, Overlap> chunkGraph = new DefaultDirectedGraph<Read, Overlap>(Overlap.class);

            Read v = startPoint;
            chunkGraph.addVertex(v);

            while (visibleOutDegreeOf(g, v) == 1) {
                for (Overlap e : g.outgoingEdgesOf(v)) {
                    Read t = g.getEdgeTarget(e);

                    chunkGraph.addVertex(t);
                    chunkGraph.addEdge(v, t, e);
                }

                Overlap e = visibleOutgoingEdgesOf(g, v).iterator().next();
                v = g.getEdgeTarget(e);
            }

            chunks.add(chunkGraph);
        }

        return chunks;
    }

    private class AlignmentGrid {
        private Map<Read, Integer> positions = new LinkedHashMap<Read, Integer>();
        private int offset = 0;

        public void addRead(Read v) {
            positions.put(v, 0);
        }

        public void addRead(Read v, Read t, AlignmentEntry ae) {
            if (v.isReversed) {
                if (ae.readA == v.getReadId()) {
                    addRead(v, t, v.getReadLength() - ae.endA);
                } else {
                    addRead(v, t, ae.endB - ae.endA);
                }
            } else {
                if (ae.readA == v.getReadId()) {
                    addRead(v, t, ae.startA - ae.startB);
                } else {
                    addRead(v, t, ae.startB - ae.startA);
                }
            }
        }

        public void addRead(Read v, Read t, int relativePosition) {
            if (!positions.containsKey(v)) {
                addRead(v);
            }

            if (!positions.containsKey(t)) {
                int sourcePos = positions.get(v);
                int absolutePosition = sourcePos + relativePosition;

                positions.put(t, absolutePosition);

                if (absolutePosition < offset) {
                    offset = absolutePosition;
                }
            }
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();

            for (Read v : positions.keySet()) {
                int position = positions.get(v);
                int absolutePosition = position - offset;

                sb.append(String.format("(%-10s) ", v.getLabel())).append(StringUtil.repeatCharNTimes(' ', absolutePosition)).append(v.getReadSequence()).append("\n");
            }

            return sb.toString();
        }
    }

    private void printChunk(DirectedGraph<Read, Overlap> chunk) {
        Read startPoint = null;

        for (Read v : chunk.vertexSet()) {
            if (visibleInDegreeOf(chunk, v) == 0 && visibleOutDegreeOf(chunk, v) == 1) {
                startPoint = v;
                break;
            }
        }

        if (startPoint != null) {
            AlignmentGrid ag = new AlignmentGrid();

            Read v = startPoint;
            while (visibleOutDegreeOf(chunk, v) == 1) {
                for (Overlap e : chunk.outgoingEdgesOf(v)) {
                    Read s = chunk.getEdgeSource(e);
                    Read t = chunk.getEdgeTarget(e);
                    AlignmentEntry ae = e.getAlignment();

                    ag.addRead(s, t, ae);

                    log.info("align ag={} s={} t={} ae={} num={}\n{}", ag.offset, s.isReversed, t.isReversed, ae, ag.positions.size(), ag);
                }

                Overlap e = visibleOutgoingEdgesOf(chunk, v).iterator().next();
                v = chunk.getEdgeTarget(e);
            }

            log.info("\n{}", ag);
        }
    }

    @Override
    public void execute() {
        log.info("Building overlap graph...");
        DirectedGraph<Read, Overlap> g = constructOverlapGraph();

        try {
            printGraph(g, new PrintStream(out), false);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        markTransitivelyInferribleEdges(g, 3);
        markAmbiguousEdgesToCut(g);

        Set<DirectedGraph<Read, Overlap>> chunks = getChunks(g);

        log.info("  stats");
        log.info("  - nodes:  {}", g.vertexSet().size());
        log.info("  - edges:  {}", g.edgeSet().size());
        log.info("  - chunks: {}", chunks.size());

        log.info("Chunks:");
        int chunkIndex = 0;
        for (DirectedGraph<Read, Overlap> chunk : chunks) {
            log.info("  chunk {}:", chunkIndex);
            chunkIndex++;

            printChunk(chunk);
            printChunk(chunk);
        }

        try {
            printGraph(g, new PrintStream(rout), true);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}
