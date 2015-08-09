package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import htsjdk.samtools.util.Interval;
import org.jgrapht.DirectedGraph;
import org.jgrapht.GraphPath;
import org.jgrapht.Graphs;
import org.jgrapht.alg.DijkstraShortestPath;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.EdgeReversedGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.containers.DataTables;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.math.MoreMathUtils;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public class GenotypeGraph extends Module {
    @Argument(fullName = "graph", shortName = "g", doc = "Graph")
    public CortexGraph GRAPH;

    @Argument(fullName = "graphRaw", shortName = "r", doc = "Graph (raw)", required = false)
    public CortexGraph GRAPH_RAW;

    @Argument(fullName = "novelGraph", shortName = "n", doc = "Graph of novel kmers")
    public CortexGraph NOVEL;

    @Argument(fullName = "ref1", shortName = "r1", doc = "Fasta file for first parent")
    public File REF1;

    @Argument(fullName = "ref2", shortName = "r2", doc = "Fasta file for second parent")
    public File REF2;

    @Argument(fullName = "bed", shortName = "b", doc = "Bed file describing variants", required = false)
    public File BED;

    @Argument(fullName = "novelKmerMap", shortName = "m", doc = "Novel kmer map", required = false)
    public File NOVEL_KMER_MAP;

    @Argument(fullName = "beAggressive", shortName = "a", doc = "Be aggressive in extending novel stretches")
    public Boolean AGGRESSIVE = false;

    @Output(fullName = "gout", shortName = "go", doc = "Graph out")
    public File gout;

    @Output
    public PrintStream out;

    @Output(fullName = "evalOut", shortName = "eout", doc = "Eval out")
    public PrintStream eout;

    private class VariantInfo {
        public String variantId;

        public String vchr;
        public int vstart;
        public int vstop;

        public String type;
        public String denovo;
        public String nahr;
        public int gcindex;
        public String ref;
        public String alt;
        public String leftFlank;
        public String rightFlank;

        public boolean found = false;
        public boolean matches = false;

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            VariantInfo that = (VariantInfo) o;

            if (vstart != that.vstart) return false;
            if (vstop != that.vstop) return false;
            if (gcindex != that.gcindex) return false;
            if (variantId != null ? !variantId.equals(that.variantId) : that.variantId != null) return false;
            if (vchr != null ? !vchr.equals(that.vchr) : that.vchr != null) return false;
            if (type != null ? !type.equals(that.type) : that.type != null) return false;
            if (denovo != null ? !denovo.equals(that.denovo) : that.denovo != null) return false;
            if (nahr != null ? !nahr.equals(that.nahr) : that.nahr != null) return false;
            if (ref != null ? !ref.equals(that.ref) : that.ref != null) return false;
            if (alt != null ? !alt.equals(that.alt) : that.alt != null) return false;
            if (leftFlank != null ? !leftFlank.equals(that.leftFlank) : that.leftFlank != null) return false;
            return !(rightFlank != null ? !rightFlank.equals(that.rightFlank) : that.rightFlank != null);
        }

        @Override
        public int hashCode() {
            int result = variantId != null ? variantId.hashCode() : 0;
            result = 31 * result + (vchr != null ? vchr.hashCode() : 0);
            result = 31 * result + vstart;
            result = 31 * result + vstop;
            result = 31 * result + (type != null ? type.hashCode() : 0);
            result = 31 * result + (denovo != null ? denovo.hashCode() : 0);
            result = 31 * result + (nahr != null ? nahr.hashCode() : 0);
            result = 31 * result + gcindex;
            result = 31 * result + (ref != null ? ref.hashCode() : 0);
            result = 31 * result + (alt != null ? alt.hashCode() : 0);
            result = 31 * result + (leftFlank != null ? leftFlank.hashCode() : 0);
            result = 31 * result + (rightFlank != null ? rightFlank.hashCode() : 0);
            return result;
        }

        @Override
        public String toString() {
            return "VariantInfo{" +
                    "variantId='" + variantId + '\'' +
                    ", vchr='" + vchr + '\'' +
                    ", vstart=" + vstart +
                    ", vstop=" + vstop +
                    ", type='" + type + '\'' +
                    ", denovo='" + denovo + '\'' +
                    ", nahr='" + nahr + '\'' +
                    ", gcindex=" + gcindex +
                    ", ref='" + ref + '\'' +
                    ", alt='" + alt + '\'' +
                    ", leftFlank='" + leftFlank + '\'' +
                    ", rightFlank='" + rightFlank + '\'' +
                    '}';
        }
    }

    private Map<String, VariantInfo> loadVariantInfoMap() {
        Map<String, VariantInfo> allVariants = new HashMap<String, VariantInfo>();

        if (BED != null) {
            TableReader bed = new TableReader(BED, "vchr", "vstart", "vstop", "info");

            for (Map<String, String> be : bed) {
                VariantInfo vi = new VariantInfo();

                vi.vchr = be.get("vchr");
                vi.vstart = Integer.valueOf(be.get("vstart"));
                vi.vstop = Integer.valueOf(be.get("vstop"));

                String[] kvpairs = be.get("info").split(";");
                for (String kvpair : kvpairs) {
                    String[] kv = kvpair.split("=");

                    if (kv[0].equals("id")) {
                        vi.variantId = kv[1];
                    }
                    if (kv[0].equals("type")) {
                        vi.type = kv[1];
                    }
                    if (kv[0].equals("denovo")) {
                        vi.denovo = kv[1];
                    }
                    if (kv[0].equals("nahr")) {
                        vi.nahr = kv[1];
                    }
                    if (kv[0].equals("gcindex")) {
                        vi.gcindex = Integer.valueOf(kv[1]);
                    }
                    if (kv[0].equals("ref")) {
                        vi.ref = kv[1];
                    }
                    if (kv[0].equals("alt")) {
                        vi.alt = kv[1];
                    }
                    if (kv[0].equals("left")) {
                        vi.leftFlank = kv[1];
                    }
                    if (kv[0].equals("right")) {
                        vi.rightFlank = kv[1];
                    }
                }

                allVariants.put(vi.variantId, vi);
            }
        }
        return allVariants;
    }

    private Map<CortexKmer, VariantInfo> loadNovelKmerMap() {
        Map<String, VariantInfo> allVariants = loadVariantInfoMap();

        Map<CortexKmer, VariantInfo> vis = new HashMap<CortexKmer, VariantInfo>();

        if (NOVEL_KMER_MAP != null) {
            TableReader tr = new TableReader(NOVEL_KMER_MAP);

            for (Map<String, String> te : tr) {
                VariantInfo vi = allVariants.get(te.get("variantId"));
                CortexKmer kmer = new CortexKmer(te.get("kmer"));

                vis.put(kmer, vi);
            }
        }

        return vis;
    }

    private String formatAttributes(Map<String, Object> attrs) {
        List<String> attrArray = new ArrayList<String>();

        for (String attr : attrs.keySet()) {
            String value = attrs.get(attr).toString();

            attrArray.add(attr + "=\"" + value + "\"");
        }

        return Joiner.on(" ").join(attrArray);
    }

    public void printGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, String prefix, boolean withText, boolean withPdf) {
        try {
            File f = new File(gout.getAbsolutePath() + "." + prefix + ".dot");
            File p = new File(gout.getAbsolutePath() + "." + prefix + ".pdf");

            PrintStream o = new PrintStream(f);

            String indent = "  ";

            o.println("digraph G {");

            if (withText) {
                o.println(indent + "node [ shape=box fontname=\"Courier New\" style=filled fillcolor=white ];");
            } else {
                o.println(indent + "node [ shape=point style=filled fillcolor=white ];");
            }

            for (AnnotatedVertex v : g.vertexSet()) {
                Map<String, Object> attrs = new TreeMap<String, Object>();

                if (!withText || g.inDegreeOf(v) == 0 || g.outDegreeOf(v) == 0) {
                    attrs.put("label", "");
                }

                if (v.isNovel()) {
                    attrs.put("color", "red");
                    attrs.put("fillcolor", "red");
                    attrs.put("shape", "circle");
                }

                if (v.flagIsSet("start") || v.flagIsSet("end")) {
                    attrs.put("label", v.getKmer());
                    attrs.put("fillcolor", v.flagIsSet("start") ? "orange" : "purple");
                    attrs.put("shape", "rect");
                }

                o.println(indent + "\"" + v.getKmer() + "\" [ " + formatAttributes(attrs) + " ];");
            }

            String[] colors = new String[]{"red", "blue", "green"};

            for (AnnotatedEdge e : g.edgeSet()) {
                String s = g.getEdgeSource(e).getKmer();
                String t = g.getEdgeTarget(e).getKmer();

                for (int c = 0; c < 3; c++) {
                    if (e.isPresent(c)) {
                        Map<String, Object> attrs = new TreeMap<String, Object>();
                        attrs.put("color", colors[c]);

                        o.println(indent + "\"" + s + "\" -> \"" + t + "\" [ " + formatAttributes(attrs) + " ];");
                    }
                }
            }

            o.println("}");

            o.close();

            if (withPdf) {
                //log.info("dot -Tpdf -o" + p.getAbsolutePath() + " " + f.getAbsolutePath());
                Runtime.getRuntime().exec("dot -Tpdf -o" + p.getAbsolutePath() + " " + f.getAbsolutePath());
            }
        } catch (FileNotFoundException e) {
            throw new IndianaException("File not found", e);
        } catch (IOException e) {
            throw new IndianaException("IO exception", e);
        }
    }

    private void addGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int color, Map<CortexKmer, Boolean> novelKmers) {
        for (AnnotatedVertex v : g.vertexSet()) {
            AnnotatedVertex av = new AnnotatedVertex(v.getKmer(), novelKmers.containsKey(new CortexKmer(v.getKmer())));

            a.addVertex(av);
        }

        for (AnnotatedEdge e : g.edgeSet()) {
            AnnotatedVertex s0 = g.getEdgeSource(e);
            AnnotatedVertex s1 = g.getEdgeTarget(e);

            CortexKmer ck0 = new CortexKmer(s0.getKmer());
            CortexKmer ck1 = new CortexKmer(s1.getKmer());

            AnnotatedVertex a0 = new AnnotatedVertex(s0.getKmer(), novelKmers.containsKey(ck0));
            AnnotatedVertex a1 = new AnnotatedVertex(s1.getKmer(), novelKmers.containsKey(ck1));

            if (!a.containsEdge(a0, a1)) {
                a.addEdge(a0, a1, new AnnotatedEdge());
            }

            a.getEdge(a0, a1).set(color, true);
        }
    }

    private DirectedGraph<AnnotatedVertex, AnnotatedEdge> copyGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a) {
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> b = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);

        for (AnnotatedVertex av : a.vertexSet()) {
            AnnotatedVertex newAv = new AnnotatedVertex(av.getKmer(), av.isNovel());

            newAv.setFlags(av.getFlags());

            b.addVertex(newAv);
        }

        for (AnnotatedEdge ae : a.edgeSet()) {
            AnnotatedVertex vs = a.getEdgeSource(ae);
            AnnotatedVertex vt = a.getEdgeTarget(ae);

            b.addEdge(vs, vt, new AnnotatedEdge(ae.getPresence()));
        }

        return b;
    }

    private DirectedGraph<AnnotatedVertex, AnnotatedEdge> simplifyGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a) {
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> b = copyGraph(a);

        AnnotatedVertex thisVertex;
        Set<AnnotatedVertex> usedStarts = new HashSet<AnnotatedVertex>();

        Set<AnnotatedVertex> candidateStarts = new HashSet<AnnotatedVertex>();
        Set<AnnotatedVertex> candidateEnds = new HashSet<AnnotatedVertex>();
        for (AnnotatedVertex av : b.vertexSet()) {
            if (av.flagIsSet("start")) {
                candidateStarts.add(av);
            } else if (av.flagIsSet("end")) {
                candidateEnds.add(av);
            }
        }

        do {
            thisVertex = null;

            Iterator<AnnotatedVertex> vertices = b.vertexSet().iterator();

            while (thisVertex == null && vertices.hasNext()) {
                AnnotatedVertex tv = vertices.next();

                if (!usedStarts.contains(tv)) {
                    while (b.inDegreeOf(tv) == 1) {
                        AnnotatedVertex pv = b.getEdgeSource(b.incomingEdgesOf(tv).iterator().next());

                        if (b.outDegreeOf(pv) == 1 && pv.isNovel() == tv.isNovel()) {
                            tv = pv;
                        } else {
                            break;
                        }
                    }

                    thisVertex = tv;
                }
            }

            if (thisVertex != null) {
                usedStarts.add(thisVertex);

                AnnotatedVertex nextVertex = b.outDegreeOf(thisVertex) == 0 ? null : b.getEdgeTarget(b.outgoingEdgesOf(thisVertex).iterator().next());

                while (nextVertex != null && b.outDegreeOf(thisVertex) == 1 && b.inDegreeOf(nextVertex) == 1 && thisVertex.isNovel() == nextVertex.isNovel()) {
                    AnnotatedVertex sv = new AnnotatedVertex(thisVertex.getKmer() + nextVertex.getKmer().charAt(nextVertex.getKmer().length() - 1), thisVertex.isNovel());

                    b.addVertex(sv);

                    Set<AnnotatedEdge> edgesToRemove = new HashSet<AnnotatedEdge>();

                    for (AnnotatedEdge e : b.incomingEdgesOf(thisVertex)) {
                        AnnotatedVertex pv = b.getEdgeSource(e);

                        b.addEdge(pv, sv, new AnnotatedEdge(e.getPresence()));

                        edgesToRemove.add(e);
                    }

                    for (AnnotatedEdge e : b.outgoingEdgesOf(nextVertex)) {
                        AnnotatedVertex nv = b.getEdgeTarget(e);

                        b.addEdge(sv, nv, new AnnotatedEdge(e.getPresence()));

                        edgesToRemove.add(e);
                    }

                    b.removeVertex(thisVertex);
                    b.removeVertex(nextVertex);
                    b.removeAllEdges(edgesToRemove);

                    thisVertex = sv;
                    nextVertex = b.outDegreeOf(thisVertex) == 0 ? null : b.getEdgeTarget(b.outgoingEdgesOf(thisVertex).iterator().next());
                }
            }
        } while (thisVertex != null);

        for (AnnotatedVertex av : b.vertexSet()) {
            for (AnnotatedVertex candidateStart : candidateStarts) {
                if (av.getKmer().contains(candidateStart.getKmer())) {
                    av.setFlag("start");
                }
            }

            for (AnnotatedVertex candidateEnd : candidateEnds) {
                if (av.getKmer().contains(candidateEnd.getKmer())) {
                    av.setFlag("end");
                }
            }
        }

        return b;
    }

    private DirectedGraph<AnnotatedVertex, AnnotatedEdge> removeOtherColors(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, int colorToRetain) {
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> b = copyGraph(a);

        Set<AnnotatedEdge> edgesToRemove = new HashSet<AnnotatedEdge>();

        for (AnnotatedEdge ae : b.edgeSet()) {
            for (int c = 0; c < 10; c++) {
                if (c != colorToRetain) {
                    ae.set(c, false);
                }
            }

            if (ae.isAbsent(colorToRetain)) {
                edgesToRemove.add(ae);
            }
        }

        b.removeAllEdges(edgesToRemove);

        Set<AnnotatedVertex> verticesToRemove = new HashSet<AnnotatedVertex>();
        for (AnnotatedVertex av : b.vertexSet()) {
            if (b.inDegreeOf(av) == 0 && b.outDegreeOf(av) == 0) {
                verticesToRemove.add(av);
            }
        }

        b.removeAllVertices(verticesToRemove);

        return b;
    }

    private String linearizePath(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, GraphPath<AnnotatedVertex, AnnotatedEdge> p) {
        StringBuilder sb = new StringBuilder(p.getStartVertex().getKmer());

        for (AnnotatedEdge ae : p.getEdgeList()) {
            AnnotatedVertex av = a.getEdgeTarget(ae);

            sb.append(av.getKmer().charAt(av.getKmer().length() - 1));
        }

        return sb.toString();
    }

    private List<String> getNovelStretch(String stretch, Map<CortexKmer, Boolean> novelKmers) {
        boolean[] novel = new boolean[stretch.length()];
        int kmerLength = novelKmers.keySet().iterator().next().length();

        for (int i = 0; i <= stretch.length() - kmerLength; i++) {
            String sk = stretch.substring(i, i + kmerLength);
            CortexKmer ck = new CortexKmer(sk);

            novel[i] = novelKmers.containsKey(ck);
        }

        List<String> novelStretches = new ArrayList<String>();

        StringBuilder novelStretchBuilder = new StringBuilder();

        for (int i = 0; i <= stretch.length() - kmerLength; i++) {
            String sk = stretch.substring(i, i + kmerLength);
            boolean isNovel = novel[i];

            if (!isNovel) {
                if (novelStretchBuilder.length() > 0) {
                    novelStretches.add(novelStretchBuilder.toString());
                    novelStretchBuilder = new StringBuilder();
                }
            } else {
                if (novelStretchBuilder.length() == 0) {
                    novelStretchBuilder.append(sk);
                } else {
                    novelStretchBuilder.append(sk.charAt(sk.length() - 1));
                }
            }
        }

        if (novelStretchBuilder.length() > 0) {
            novelStretches.add(novelStretchBuilder.toString());
        }

        return novelStretches;
    }

    private class PathInfo {
        public String start, stop, child, parent;

        public PathInfo(String start, String stop, String child, String parent) {
            this.start = start;
            this.stop = stop;
            this.child = child;
            this.parent = parent;
        }
    }

    private PathInfo computeBestMinWeightPath(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, int color, String stretch, Map<CortexKmer, Boolean> novelKmers) {
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> b = copyGraph(a);

        AnnotateStartsAndEnds annotateStartsAndEnds = new AnnotateStartsAndEnds(color, stretch, novelKmers, b).invoke();
        Set<AnnotatedVertex> candidateStarts = annotateStartsAndEnds.getCandidateStarts();
        Set<AnnotatedVertex> candidateEnds = annotateStartsAndEnds.getCandidateEnds();

        DirectedGraph<AnnotatedVertex, AnnotatedEdge> b0 = removeOtherColors(b, 0);
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> bc = removeOtherColors(b, color);

        List<String> novelStretches = getNovelStretch(stretch, novelKmers);

        double minPl0 = Double.MAX_VALUE, minPlc = Double.MAX_VALUE;
        String minLp0 = "", minLpc = "";
        String start = "", stop = "";

        for (AnnotatedVertex sv : candidateStarts) {
            for (AnnotatedVertex ev : candidateEnds) {
                DijkstraShortestPath<AnnotatedVertex, AnnotatedEdge> dsp0 = new DijkstraShortestPath<AnnotatedVertex, AnnotatedEdge>(b0, sv, ev);
                DijkstraShortestPath<AnnotatedVertex, AnnotatedEdge> dspc = new DijkstraShortestPath<AnnotatedVertex, AnnotatedEdge>(bc, sv, ev);

                GraphPath<AnnotatedVertex, AnnotatedEdge> p0 = dsp0.getPath();
                GraphPath<AnnotatedVertex, AnnotatedEdge> pc = dspc.getPath();

                String lp0 = p0 == null ? "" : linearizePath(b0, p0);
                String lpc = pc == null ? "" : linearizePath(bc, pc);

                if (lp0.contains(sv.getKmer()) && lp0.contains(ev.getKmer()) && lpc.contains(sv.getKmer()) && lpc.contains(ev.getKmer())) {
                    boolean stretchesArePresent = true;
                    for (String novelStretch : novelStretches) {
                        if (!lp0.contains(novelStretch)) {
                            stretchesArePresent = false;
                            break;
                        }
                    }

                    if (novelStretches.size() > 0 && stretchesArePresent) {
                        double pl0 = dsp0.getPathLength();
                        double plc = dspc.getPathLength();

                        if (pl0 < minPl0 && plc < minPlc) {
                            minPl0 = pl0;
                            minLp0 = lp0;

                            minPlc = plc;
                            minLpc = lpc;

                            start = sv.getKmer();
                            stop = ev.getKmer();
                        }
                    }
                }
            }
        }

        return new PathInfo(start, stop, minLp0, minLpc);
    }

    private boolean hasRecombinations(String stretch) {
        StringBuilder inherit = new StringBuilder();
        for (int i = 0; i <= stretch.length() - GRAPH.getKmerSize(); i++) {
            String sk = stretch.substring(i, i + GRAPH.getKmerSize());
            CortexKmer ck = new CortexKmer(sk);

            CortexRecord cr = GRAPH.findRecord(ck);

            if (cr != null) {
                int cov0 = cr.getCoverage(0);
                int cov1 = cr.getCoverage(1);
                int cov2 = cr.getCoverage(2);

                if (cov0 > 0 && cov1 == 0 && cov2 == 0) {
                    inherit.append("C");
                } else if (cov0 > 0 && cov1 > 0 && cov2 == 0) {
                    inherit.append("M");
                } else if (cov0 > 0 && cov1 == 0 && cov2 > 0) {
                    inherit.append("D");
                } else if (cov0 > 0 && cov1 > 0 && cov2 > 0) {
                    inherit.append("B");
                } else {
                    inherit.append("?");
                }
            }
        }

        for (int i = 0; i < inherit.length(); i++) {
            if (inherit.charAt(i) == 'B') {
                char prevContext = '?';
                char nextContext = '?';

                for (int j = i - 1; j >= 0; j--) {
                    if (inherit.charAt(j) == 'M' || inherit.charAt(j) == 'D') {
                        prevContext = inherit.charAt(j);
                        break;
                    }
                }

                for (int j = i + 1; j < inherit.length(); j++) {
                    if (inherit.charAt(j) == 'M' || inherit.charAt(j) == 'D') {
                        nextContext = inherit.charAt(j);
                        break;
                    }
                }

                char context = '?';
                if (prevContext == nextContext && prevContext != '?') {
                    context = prevContext;
                } else if (prevContext != nextContext) {
                    if (prevContext != '?') {
                        context = prevContext;
                    } else {
                        context = nextContext;
                    }
                }

                inherit.setCharAt(i, context);
            }
        }

        log.debug("    - inherit:     {}", inherit.toString());

        String inheritStr = inherit.toString();
        return inheritStr.matches(".*D+.C+.M+.*") || inheritStr.matches(".*M+.C+.D+.*");
    }

    private boolean isChimeric(String stretch, KmerLookup kl) {
        List<Set<Interval>> ils = kl.findKmers(stretch);

        StringBuilder sb = new StringBuilder();

        int nextAvailableCode = 0;
        Map<String, Integer> chrCodes = new HashMap<String, Integer>();
        Map<String, Integer> chrCount = new HashMap<String, Integer>();

        for (int i = 0; i < ils.size(); i++) {
            Set<Interval> il = ils.get(i);

            if (il.size() == 1) {
                Interval interval = il.iterator().next();

                ContainerUtils.increment(chrCount, interval.getSequence());

                if (!chrCodes.containsKey(interval.getSequence())) {
                    chrCodes.put(interval.getSequence(), nextAvailableCode);
                    nextAvailableCode++;
                }

                int code = chrCodes.get(interval.getSequence());
                sb.append(code);
            } else {
                sb.append(".");
            }
        }

        log.debug("    - chimeric:    {}", sb.toString());

        return chrCount.size() > 1;
    }

    private GraphicalVariantContext callVariant(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, int color, String stretch, Map<CortexKmer, Boolean> novelKmers, KmerLookup kl) {
        // Compute paths
        PathInfo p = computeBestMinWeightPath(a, color, stretch, novelKmers);

        // Trim back to reference and variant alleles
        int s, e0 = p.child.length() - 1, e1 = p.parent.length() - 1;

        for (s = 0; s < (p.child.length() < p.parent.length() ? p.child.length() : p.parent.length()) && p.child.charAt(s) == p.parent.charAt(s); s++) {
        }

        while (e0 > s && e1 > s && p.child.charAt(e0) == p.parent.charAt(e1)) {
            e0--;
            e1--;
        }

        String parentalAllele = p.parent == null || p.parent.equals("") || p.child.equals(p.parent) ? "A" : p.parent.substring(s, e1 + 1);
        String childAllele = p.child == null || p.child.equals("") || p.child.equals(p.parent) ? "N" : p.child.substring(s, e0 + 1);

        int e = s + parentalAllele.length() - 1;

        // Decide if the event is actually a GC or NAHR event
        log.debug("    novel stretch: {}", stretch);

        boolean hasRecombs = hasRecombinations(stretch);
        boolean isChimeric = isChimeric(stretch, kl);

        List<Set<Interval>> alignment = kl.align(CortexUtils.getSeededStretchLeft(GRAPH, p.start, color, false) + p.parent + CortexUtils.getSeededStretchRight(GRAPH, p.stop, color, false));
        List<Set<Interval>> anovel = kl.align(stretch);

        // Build the GVC
        GraphicalVariantContext gvc = new GraphicalVariantContext()
                .attribute(color, "start", s)
                .attribute(color, "stop", e)
                .attribute(color, "parentalAllele", parentalAllele)
                .attribute(color, "childAllele", childAllele)
                .attribute(color, "parentalStretch", p.parent)
                .attribute(color, "childStretch", p.child)
                .attribute(color, "event", "unknown")
                .attribute(color, "traversalStatus", "complete")
                .attribute(color, "parentalPathAlignment", alignment)
                .attribute(color, "novelStretchAlignment", anovel)
                .attribute(color, "haplotypeBackground", color);

        if (childAllele.equals("N")) {
            if (hasRecombs) {
                gvc.attribute(color, "event", "GC");
            } else if (isChimeric) {
                gvc.attribute(color, "event", "NAHR");
            } else {
                gvc.attribute(color, "traversalStatus", "incomplete");
            }
        } else {
            if (parentalAllele.length() == childAllele.length()) {
                if (parentalAllele.length() == 1) {
                    gvc.attribute(color, "event", "SNP");
                } else if (SequenceUtils.reverseComplement(parentalAllele).equals(childAllele)) {
                    gvc.attribute(color, "event", "INV");
                }
            } else if (parentalAllele.length() == 1 && childAllele.length() > 1) {
                gvc.attribute(color, "event", "INS");

                String childAlleleTrimmed = childAllele;
                if (parentalAllele.length() == 1 && childAllele.charAt(0) == parentalAllele.charAt(0)) {
                    childAlleleTrimmed = childAllele.substring(1, childAllele.length());
                } else if (parentalAllele.length() == 1 && childAllele.charAt(childAllele.length() - 1) == parentalAllele.charAt(0)) {
                    childAlleleTrimmed = childAllele.substring(0, childAllele.length() - 1);
                }

                String repUnit = "";
                boolean isStrExp = false;
                for (int repLength = 2; repLength <= 5 && repLength <= childAlleleTrimmed.length(); repLength++) {
                    repUnit = childAlleleTrimmed.substring(0, repLength);

                    boolean fits = true;
                    for (int i = 0; i < childAlleleTrimmed.length(); i += repLength) {
                        fits &= (i + repLength <= childAllele.length()) && childAllele.substring(i, i + repLength).equals(repUnit);
                    }

                    if (fits &&
                            (p.parent.substring(s - repLength, s).equals(repUnit) ||
                                    p.parent.substring(s, s + repLength).equals(repUnit))) {
                        gvc.attribute(color, "event", "STR_EXP");
                        gvc.attribute(color, "repeatingUnit", repUnit);
                        isStrExp = true;
                    }
                }

                if (!isStrExp && childAlleleTrimmed.length() >= 10 && childAlleleTrimmed.length() <= 50 &&
                        ((s - childAlleleTrimmed.length() >= 0 && p.parent.substring(s - childAlleleTrimmed.length(), s).equals(childAlleleTrimmed)) ||
                                (s + childAlleleTrimmed.length() <= p.parent.length() && p.parent.substring(s, s + childAlleleTrimmed.length()).equals(childAlleleTrimmed)))) {
                    gvc.attribute(color, "event", "TD");
                    gvc.attribute(color, "repeatingUnit", repUnit);
                }
            } else if (parentalAllele.length() > 1 && childAllele.length() == 1) {
                gvc.attribute(color, "event", "DEL");

                String parentalAlleleTrimmed = parentalAllele;
                if (childAllele.length() == 1 && parentalAllele.charAt(0) == childAllele.charAt(0)) {
                    parentalAlleleTrimmed = parentalAllele.substring(1, parentalAllele.length());
                } else if (childAllele.length() == 1 && parentalAllele.charAt(parentalAllele.length() - 1) == childAllele.charAt(0)) {
                    parentalAlleleTrimmed = parentalAllele.substring(0, parentalAllele.length() - 1);
                }

                String repUnit = "";
                for (int repLength = 2; repLength <= 5 && repLength <= parentalAlleleTrimmed.length(); repLength++) {
                    repUnit = parentalAlleleTrimmed.substring(0, repLength);

                    boolean fits = true;
                    for (int i = 0; i < parentalAlleleTrimmed.length(); i += repLength) {
                        fits &= (i + repLength <= parentalAllele.length()) && parentalAllele.substring(i, i + repLength).equals(repUnit);
                    }

                    if (fits &&
                            (p.child.substring(s - repLength, s).equals(repUnit) ||
                                    p.child.substring(s + parentalAlleleTrimmed.length(), s + parentalAlleleTrimmed.length() + repLength).equals(repUnit))) {
                        gvc.attribute(color, "event", "STR_CON");
                        gvc.attribute(color, "repeatingUnit", repUnit);
                    }
                }
            }
        }

        return gvc;
    }

    private DirectedGraph<AnnotatedVertex, AnnotatedEdge> loadLocalGraph(Map<CortexKmer, Boolean> novelKmers, String stretch) {
        // first, explore each color and bring the local subgraphs into memory
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg0 = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg1 = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg2 = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);

        for (int i = 0; i <= stretch.length() - GRAPH.getKmerSize(); i++) {
            String kmer = stretch.substring(i, i + GRAPH.getKmerSize());
            AnnotatedVertex ak = new AnnotatedVertex(kmer);

            if (!sg0.containsVertex(ak)) {
                Graphs.addGraph(sg0, CortexUtils.dfs(GRAPH, GRAPH_RAW, kmer, 0, null, new AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge>() {
                    @Override
                    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int depth) {
                        return cr.getCoverage(1) > 0 || cr.getCoverage(2) > 0;
                    }

                    @Override
                    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions) {
                        return false;
                    }

                    @Override
                    public int maxJunctionsAllowed() {
                        return 3;
                    }
                }));
            }
        }

        for (int c = 1; c <= 2; c++) {
            for (AnnotatedVertex ak : sg0.vertexSet()) {
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg = (c == 1) ? sg1 : sg2;

                if (!sg.containsVertex(ak)) {
                    TraversalStopper<AnnotatedVertex, AnnotatedEdge> stopper = new AbstractTraversalStopper<AnnotatedVertex, AnnotatedEdge>() {
                        @Override
                        public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions) {
                            String fw = cr.getKmerAsString();
                            String rc = SequenceUtils.reverseComplement(fw);

                            if (g.containsVertex(new AnnotatedVertex(fw)) || g.containsVertex(new AnnotatedVertex(rc))) {
                                if (junctions < distanceToGoal) {
                                    distanceToGoal = junctions;
                                }

                                return true;
                            }

                            return false;
                        }

                        @Override
                        public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int junctions) {
                            return junctions >= maxJunctionsAllowed();
                        }

                        @Override
                        public int maxJunctionsAllowed() {
                            return 5;
                        }
                    };

                    Graphs.addGraph(sg, CortexUtils.dfs(GRAPH, GRAPH_RAW, ak.getKmer(), c, sg0, stopper));
                }
            }
        }

        // Now, combine them all into an annotated graph
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> ag = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);

        addGraph(ag, sg0, 0, novelKmers);
        addGraph(ag, sg1, 1, novelKmers);
        addGraph(ag, sg2, 2, novelKmers);

        return ag;
    }

    private class AnnotateStartsAndEnds {
        private int color;
        private String stretch;
        private Map<CortexKmer, Boolean> novelKmers;
        private DirectedGraph<AnnotatedVertex, AnnotatedEdge> b;
        private Set<AnnotatedVertex> candidateEnds;
        private Set<AnnotatedVertex> candidateStarts;

        public AnnotateStartsAndEnds(int color, String stretch, Map<CortexKmer, Boolean> novelKmers, DirectedGraph<AnnotatedVertex, AnnotatedEdge> b) {
            this.color = color;
            this.stretch = stretch;
            this.novelKmers = novelKmers;
            this.b = b;
        }

        public Set<AnnotatedVertex> getCandidateEnds() {
            return candidateEnds;
        }

        public Set<AnnotatedVertex> getCandidateStarts() {
            return candidateStarts;
        }

        public AnnotateStartsAndEnds invoke() {
            int kmerLength = novelKmers.keySet().iterator().next().length();

            String fk = stretch.substring(0, kmerLength);
            AnnotatedVertex afk = new AnnotatedVertex(fk, novelKmers.containsKey(new CortexKmer(fk)));

            String lk = stretch.substring(stretch.length() - kmerLength, stretch.length());
            AnnotatedVertex alk = new AnnotatedVertex(lk, novelKmers.containsKey(new CortexKmer(lk)));

            candidateEnds = new HashSet<AnnotatedVertex>();
            DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> dfsLk = new DepthFirstIterator<AnnotatedVertex, AnnotatedEdge>(b, afk);
            while (dfsLk.hasNext()) {
                AnnotatedVertex av = dfsLk.next();

                if (b.inDegreeOf(av) > 1) {
                    Set<AnnotatedEdge> aes = b.incomingEdgesOf(av);

                    boolean adjacentToNovelty = false;
                    AnnotatedEdge edgeToNovelty = null;

                    for (AnnotatedEdge ae : aes) {
                        if (ae.isPresent(0) && ae.isAbsent(1) && ae.isAbsent(2)) {
                            adjacentToNovelty = true;
                            edgeToNovelty = ae;
                        }
                    }

                    if (adjacentToNovelty) {
                        for (AnnotatedEdge ae : aes) {
                            if (!ae.equals(edgeToNovelty)) {
                                if (ae.isPresent(color)) {
                                    candidateEnds.add(av);

                                    //av.setFlag("end");

                                    //
                                    for (AnnotatedVertex ava : b.vertexSet()) {
                                        if (av.equals(ava)) {
                                            ava.setFlag("end");
                                        }
                                    }
                                    //
                                }
                            }
                        }
                    }
                }
            }

            candidateStarts = new HashSet<AnnotatedVertex>();
            DirectedGraph<AnnotatedVertex, AnnotatedEdge> rag = new EdgeReversedGraph<AnnotatedVertex, AnnotatedEdge>(b);
            DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> dfsFk = new DepthFirstIterator<AnnotatedVertex, AnnotatedEdge>(rag, alk);
            while (dfsFk.hasNext()) {
                AnnotatedVertex av = dfsFk.next();

                if (rag.inDegreeOf(av) > 1) {
                    Set<AnnotatedEdge> aes = rag.incomingEdgesOf(av);

                    boolean adjacentToNovelty = false;
                    AnnotatedEdge edgeToNovelty = null;

                    for (AnnotatedEdge ae : aes) {
                        if (ae.isPresent(0) && ae.isAbsent(1) && ae.isAbsent(2)) {
                            adjacentToNovelty = true;
                            edgeToNovelty = ae;
                        }
                    }

                    if (adjacentToNovelty) {
                        for (AnnotatedEdge ae : aes) {
                            if (!ae.equals(edgeToNovelty)) {
                                if (ae.isPresent(color)) {
                                    candidateStarts.add(av);

                                    //av.setFlag("start");

                                    //
                                    for (AnnotatedVertex ava : b.vertexSet()) {
                                        if (av.equals(ava)) {
                                            ava.setFlag("start");
                                        }
                                    }
                                    //
                                }
                            }
                        }
                    }
                }
            }
            return this;
        }
    }

    private void evalVariant(GraphicalVariantContext gvc, int color, Map<CortexKmer, VariantInfo> vis, String stretch) {
        Set<VariantInfo> relevantVis = new HashSet<VariantInfo>();
        int kmerSize = vis.keySet().iterator().next().length();

        for (int i = 0; i <= stretch.length() - kmerSize; i++) {
            CortexKmer ck = new CortexKmer(stretch.substring(i, i + kmerSize));

            if (vis.containsKey(ck)) {
                relevantVis.add(vis.get(ck));
            }
        }

        if (relevantVis.size() > 0) {
            VariantInfo bestVi = null;

            for (VariantInfo vi : relevantVis) {
                String ref = gvc.getAttributeAsString(color, "parentalAllele");
                String alt = gvc.getAttributeAsString(color, "childAllele");

                String refStretch = gvc.getAttributeAsString(color, "parentalStretch");
                String altStretch = gvc.getAttributeAsString(color, "childStretch");

                int pos = gvc.getAttributeAsInt(color, "start");
                int refLength = gvc.getAttributeAsString(color, "parentalAllele").length();
                int altLength = gvc.getAttributeAsString(color, "childAllele").length();
                boolean found = false;

                while (pos >= 0 && pos + refLength < refStretch.length() && pos + altLength < altStretch.length()) {
                    String refFw = refStretch.substring(pos, pos + refLength);
                    String refRc = SequenceUtils.reverseComplement(refFw);

                    String altFw = altStretch.substring(pos, pos + altLength);
                    String altRc = SequenceUtils.reverseComplement(altFw);

                    if (vi.ref.equals(refFw) && vi.alt != null && vi.alt.equals(altFw)) {
                        ref = refFw;
                        alt = altFw;
                        found = true;
                        break;
                    } else if (vi.ref.equals(refRc) && vi.alt != null && vi.alt.equals(altRc)) {
                        ref = refRc;
                        alt = altRc;
                        found = true;
                        break;
                    }

                    pos--;
                }

                if (!found) {
                    pos = gvc.getAttributeAsInt(color, "start");

                    while (pos >= 0 && pos + refLength < refStretch.length() && pos + altLength < altStretch.length()) {
                        String refFw = refStretch.substring(pos, pos + refLength);
                        String refRc = SequenceUtils.reverseComplement(refFw);

                        String altFw = altStretch.substring(pos, pos + altLength);
                        String altRc = SequenceUtils.reverseComplement(altFw);

                        if (vi.ref.equals(refFw) && vi.alt != null && vi.alt.equals(altFw)) {
                            ref = refFw;
                            alt = altFw;
                            found = true;
                            break;
                        } else if (vi.ref.equals(refRc) && vi.alt != null && vi.alt.equals(altRc)) {
                            ref = refRc;
                            alt = altRc;
                            found = true;
                            break;
                        }

                        pos++;
                    }
                }

                String knownRef = vi.ref;
                String knownAlt = vi.alt == null ? "" : vi.alt;

                if (!found && knownRef.length() > ref.length() && knownAlt.length() > alt.length() && knownRef.length() == knownAlt.length()) {
                    String refFw = ref;
                    String altFw = alt;
                    String refRc = SequenceUtils.reverseComplement(refFw);
                    String altRc = SequenceUtils.reverseComplement(altFw);

                    String refFinal = null, altFinal = null;

                    if (knownRef.contains(refFw) && knownAlt.contains(altFw)) {
                        refFinal = refFw;
                        altFinal = altFw;
                    } else if (knownRef.contains(refRc) && knownAlt.contains(altRc)) {
                        refFinal = refRc;
                        altFinal = altRc;
                    }

                    if (refFinal != null && altFinal != null) {
                        int r0index = knownRef.indexOf(refFinal);
                        int a0index = knownAlt.indexOf(altFinal);
                        int r1index = r0index + refFinal.length();
                        int a1index = a0index + altFinal.length();

                        if (r0index == a0index && r1index == a1index &&
                                knownRef.substring(0, r0index).equals(knownAlt.substring(0, a0index)) &&
                                knownRef.substring(r1index, knownRef.length()).equals(knownAlt.substring(a1index, knownAlt.length()))) {
                            knownRef = ref;
                            knownAlt = alt;
                        }
                    }
                }

                log.debug("    - vi: {}", vi);
                log.debug("        known ref: {}", knownRef);
                log.debug("        known alt: {}", knownAlt);
                log.debug("        called ref: {}", ref);
                log.debug("        called alt: {}", alt);

                log.debug("    - matches: {} ({} {} {} {})", knownRef.equals(ref) && knownAlt.equals(alt), ref, alt, knownRef, knownAlt);

                vi.found = true;
                vi.matches = (knownRef.equals(ref) && knownAlt.equals(alt)) || vi.denovo.equals(gvc.getAttributeAsString(color, "event"));

                if (bestVi == null || vi.matches) {
                    bestVi = vi;
                }
            }

            if (bestVi != null) {
                String knownRef = bestVi.ref;
                String knownAlt = bestVi.alt;
                String ref = gvc.getAttributeAsString(color, "parentalAllele");
                String alt = gvc.getAttributeAsString(color, "childAllele");

                gvc.attribute(color, "isKnownVariant", true);
                gvc.attribute(color, "variantId", bestVi.variantId);
                gvc.attribute(color, "allelesMatch", (knownRef.equals(ref) && knownAlt.equals(alt)) || bestVi.denovo.equals(gvc.getAttributeAsString(color, "event")));
                gvc.attribute(color, "eventsMatch", bestVi.denovo.equals(gvc.getAttributeAsString(color, "event")));
            } else {
                gvc.attribute(color, "isKnownVariant", false);
                gvc.attribute(color, "variantId", "unknown");
                gvc.attribute(color, "allelesMatch", false);
                gvc.attribute(color, "eventsMatch", false);
            }
        }
    }

    public void chooseVariant(GraphicalVariantContext gvc) {
        int[] scores = new int[2];

        for (int c = 1; c <= 2; c++) {
            if (gvc.getAttribute(c, "traversalStatus").equals("complete")) { scores[c-1]++; }
            if (!gvc.getAttribute(c, "event").equals("unknown")) { scores[c-1]++; }
            if (gvc.getAttribute(c, "event").equals("GC") || gvc.getAttribute(c, "event").equals("NAHR")) { scores[c-1]++; }

            gvc.attribute(c, "score", scores[c-1]);
        }

        int bc = MoreMathUtils.whichMax(scores) + 1;

        gvc.addAttributes(0, gvc.getAttributes(bc));
        gvc.attribute(0, "haplotypeBackground", scores[0] == scores[1] ? 0 : bc);
    }

    @Override
    public void execute() {
        log.info("Loading reference indices for fast kmer lookup...");
        KmerLookup kl1 = new KmerLookup(REF1);
        KmerLookup kl2 = new KmerLookup(REF2);

        Map<CortexKmer, VariantInfo> vis = loadNovelKmerMap();
        Map<String, VariantInfo> vids = new HashMap<String, VariantInfo>();
        Map<String, Boolean> viSeen = new HashMap<String, Boolean>();
        for (VariantInfo vi : vis.values()) {
            vids.put(vi.variantId, vi);
            viSeen.put(vi.variantId, false);
        }

        Map<CortexKmer, Boolean> novelKmers = new HashMap<CortexKmer, Boolean>();
        for (CortexRecord cr : NOVEL) {
            novelKmers.put(new CortexKmer(cr.getKmerAsString()), true);
        }

        int totalNovelKmersUsed = 0;
        int stretchNum = 1;

        Set<GraphicalVariantContext> gvcs = new LinkedHashSet<GraphicalVariantContext>();
        /*
        log.info("Discovering novel stretches in graph...");
        for (CortexKmer novelKmer : novelKmers.keySet()) {
            if (novelKmers.get(novelKmer)) {
                int novelKmersUsed = 0;

                // Walk the graph left and right of novelKmer and extract a novel stretch
                String stretch = CortexUtils.getSeededStretch(GRAPH, novelKmer.getKmerAsString(), 0, AGGRESSIVE);

                // See how many novel kmers we've used up
                for (int i = 0; i <= stretch.length() - GRAPH.getKmerSize(); i++) {
                    CortexKmer ck = new CortexKmer(stretch.substring(i, i + GRAPH.getKmerSize()));

                    if (novelKmers.containsKey(ck) && novelKmers.get(ck)) {
                        totalNovelKmersUsed++;
                        novelKmersUsed++;
                        novelKmers.put(ck, false);
                    }
                }

                GraphicalVariantContext gvc = new GraphicalVariantContext()
                        .attribute(0, "stretch", stretch)
                        .attribute(0, "stretchNum", stretchNum)
                        .attribute(0, "stretchLength", stretch.length())
                        .attribute(0, "novelKmersUsed", novelKmersUsed)
                        .attribute(0, "novelKmersTotal", novelKmers.size());

                gvcs.add(gvc);

                log.info("  stretch {}: {} bp, {}/{} novel kmers, {}/{} cumulative novel kmers seen", stretchNum, stretch.length(), novelKmersUsed, novelKmers.size(), totalNovelKmersUsed, novelKmers.size());

                stretchNum++;
            }
        }
        log.info("  found {} stretches", gvcs.size());
        */

        DataTables evalTables = new DataTables();

        evalTables.addTable("variantStats", "Statistics on variants", "knownVariantId", "knownVariantEvent", "knownVariantLength", "variantId", "variantEvent", "variantLength", "novelKmersUsed");

        evalTables.addTable("discoveryStats", "Statistics on variant discovery", "tp", "tn", "fp", "fn");
        evalTables.getTable("discoveryStats").set("dummy", "tp", 0l);
        evalTables.getTable("discoveryStats").set("dummy", "tn", 0l);
        evalTables.getTable("discoveryStats").set("dummy", "fp", 0l);
        evalTables.getTable("discoveryStats").set("dummy", "fn", 0l);

        evalTables.addTable("alleleMatchStats", "Statistics on allele matches", "match", "mismatch");
        evalTables.getTable("alleleMatchStats").set("dummy", "match", 0l);
        evalTables.getTable("alleleMatchStats").set("dummy", "mismatch", 0l);

        evalTables.addTable("eventMatchStats", "Statistics on event matches", "match", "mismatch");
        evalTables.getTable("eventMatchStats").set("dummy", "match", 0l);
        evalTables.getTable("eventMatchStats").set("dummy", "mismatch", 0l);

        log.info("Genotyping novel kmer stretches in graph...");
        for (CortexKmer novelKmer : novelKmers.keySet()) {
            if (novelKmers.get(novelKmer)) {
                // Walk the graph left and right of novelKmer and extract a novel stretch
                String stretch = CortexUtils.getSeededStretch(GRAPH, novelKmer.getKmerAsString(), 0, AGGRESSIVE);

                // Construct GVC
                GraphicalVariantContext gvc = new GraphicalVariantContext()
                        .attribute(0, "stretch", stretch)
                        .attribute(0, "stretchNum", stretchNum)
                        .attribute(0, "stretchLength", stretch.length())
                        .attribute(0, "novelKmersTotal", novelKmers.size());

                log.info("  stretch {}: {} bp", stretchNum, stretch.length());

                // Fetch the local subgraph context from disk
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> ag = loadLocalGraph(novelKmers, stretch);
                log.info("    subgraph : {} vertices, {} edges", ag.vertexSet().size(), ag.edgeSet().size());

                // Extract parental stretches
                PathInfo p1 = computeBestMinWeightPath(ag, 1, stretch, novelKmers);
                PathInfo p2 = computeBestMinWeightPath(ag, 2, stretch, novelKmers);

                log.info("    paths:");
                log.info("    - 1: {}", SequenceUtils.truncate(p1.parent, 100));
                log.info("      c: {}", SequenceUtils.truncate(p1.child, 100));
                log.info("    - 2: {}", SequenceUtils.truncate(p2.parent, 100));
                log.info("      c: {}", SequenceUtils.truncate(p2.child, 100));

                // Call variants
                gvc.add(callVariant(ag, 1, stretch, novelKmers, kl1));
                gvc.add(callVariant(ag, 2, stretch, novelKmers, kl2));

                log.info("    variants:");
                log.info("    - 1: {} {} ({} bp)", gvc.getAttributeAsString(1, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(1, "parentalAllele"), 70), gvc.getAttributeAsString(1, "parentalAllele").length());
                log.info("      c: {} {} ({} bp)", gvc.getAttributeAsString(1, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(1, "childAllele"), 70), gvc.getAttributeAsString(1, "childAllele").length());
                log.info("    - 2: {} {} ({} bp)", gvc.getAttributeAsString(2, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(2, "parentalAllele"), 70), gvc.getAttributeAsString(2, "parentalAllele").length());
                log.info("      c: {} {} ({} bp)", gvc.getAttributeAsString(2, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(2, "childAllele"), 70), gvc.getAttributeAsString(2, "childAllele").length());

                // Finalize into a single call
                chooseVariant(gvc);

                // Show alignment
                log.info("    alignment:");
                log.info("    - novel stretch: {}", gvc.getAttribute(0, "novelStretchAlignment"));
                log.info("    - parental path: {}", gvc.getAttribute(0, "parentalPathAlignment"));

                // See how many novel kmers we've used up
                int novelKmersUsed = 0;

                for (int i = 0; i <= p1.child.length() - GRAPH.getKmerSize(); i++) {
                    CortexKmer ck = new CortexKmer(p1.child.substring(i, i + GRAPH.getKmerSize()));

                    if (novelKmers.containsKey(ck) && novelKmers.get(ck)) {
                        totalNovelKmersUsed++;
                        novelKmersUsed++;
                        novelKmers.put(ck, false);
                    }
                }

                for (int i = 0; i <= p2.child.length() - GRAPH.getKmerSize(); i++) {
                    CortexKmer ck = new CortexKmer(p2.child.substring(i, i + GRAPH.getKmerSize()));

                    if (novelKmers.containsKey(ck) && novelKmers.get(ck)) {
                        totalNovelKmersUsed++;
                        novelKmersUsed++;
                        novelKmers.put(ck, false);
                    }
                }

                for (int i = 0; i <= stretch.length() - GRAPH.getKmerSize(); i++) {
                    CortexKmer ck = new CortexKmer(stretch.substring(i, i + GRAPH.getKmerSize()));

                    if (novelKmers.containsKey(ck) && novelKmers.get(ck)) {
                        totalNovelKmersUsed++;
                        novelKmersUsed++;
                        novelKmers.put(ck, false);
                    }
                }

                gvc.attribute(0, "novelKmersUsed", novelKmersUsed);

                log.info("    novelty:");
                log.info("    - novel kmers used: {}/{}", novelKmersUsed, novelKmers.size());
                log.info("    - cumulative usage: {}/{}", totalNovelKmersUsed, novelKmers.size());

                // Evaluate variants
                if (BED != null) {
                    log.info("    evaluate:");

                    for (int c = 0; c < 3; c++) {
                        evalVariant(gvc, c, vis, stretch);

                        String vid = gvc.getAttributeAsString(c, "variantId");
                        VariantInfo vi = vids.containsKey(vid) ? vids.get(vid) : null;
                        String event = vi == null ? "none" : vi.denovo;

                        log.info("    - {}: background {}, isKnownVariant {}, allelesMatch {}, eventsMatch {}, variantId {} {}",
                                c,
                                gvc.getAttributeAsInt(c, "haplotypeBackground"),
                                gvc.getAttribute(c, "isKnownVariant"),
                                gvc.getAttribute(c, "allelesMatch"),
                                gvc.getAttribute(c, "eventsMatch"),
                                gvc.getAttributeAsString(c, "variantId"),
                                event
                        );
                    }

                    if (gvc.getAttributeAsBoolean(0, "isKnownVariant")) {
                        evalTables.getTable("discoveryStats").increment("dummy", "tp");

                        String vid = gvc.getAttributeAsString(0, "variantId");
                        VariantInfo vi = vids.get(vid);

                        int refLength = vi.ref != null ? vi.ref.length() : 0;
                        int altLength = vi.alt != null ? vi.alt.length() : 0;

                        String pk = vi.variantId + "." + gvc.getAttributeAsInt(0, "stretchNum");

                        evalTables.getTable("variantStats").set(pk, "knownVariantId", gvc.getAttributeAsString(0, "variantId"));
                        evalTables.getTable("variantStats").set(pk, "knownVariantEvent", vi.denovo);
                        evalTables.getTable("variantStats").set(pk, "knownVariantLength", Math.abs(refLength - altLength));
                        evalTables.getTable("variantStats").set(pk, "variantId", gvc.getAttributeAsInt(0, "stretchNum"));
                        evalTables.getTable("variantStats").set(pk, "variantEvent", gvc.getAttributeAsString(0, "event"));
                        evalTables.getTable("variantStats").set(pk, "variantLength", Math.abs(gvc.getAttributeAsString(0, "parentalAllele").length() - gvc.getAttributeAsString(0, "childAllele").length()));
                        evalTables.getTable("variantStats").set(pk, "novelKmersUsed", novelKmersUsed);

                        viSeen.put(gvc.getAttributeAsString(0, "knownVariantId"), true);

                        if (gvc.getAttributeAsBoolean(0, "allelesMatch")) {
                            evalTables.getTable("alleleMatchStats").increment("dummy", "match");
                        } else {
                            evalTables.getTable("alleleMatchStats").increment("dummy", "mismatch");
                        }

                        if (gvc.getAttributeAsBoolean(0, "eventsMatch")) {
                            evalTables.getTable("eventMatchStats").increment("dummy", "match");
                        } else {
                            evalTables.getTable("eventMatchStats").increment("dummy", "mismatch");
                        }
                    } else {
                        evalTables.getTable("discoveryStats").increment("dummy", "fp");

                        String pk = "none." + gvc.getAttributeAsInt(0, "stretchNum");

                        evalTables.getTable("variantStats").set(pk, "knownVariantId", "none");
                        evalTables.getTable("variantStats").set(pk, "knownVariantEvent", "none");
                        evalTables.getTable("variantStats").set(pk, "knownVariantLength", 0);
                        evalTables.getTable("variantStats").set(pk, "variantId", gvc.getAttributeAsInt(0, "stretchNum"));
                        evalTables.getTable("variantStats").set(pk, "variantEvent", gvc.getAttributeAsString(0, "event"));
                        evalTables.getTable("variantStats").set(pk, "variantLength", Math.abs(gvc.getAttributeAsString(0, "parentalAllele").length() - gvc.getAttributeAsString(0, "childAllele").length()));
                        evalTables.getTable("variantStats").set(pk, "novelKmersUsed", novelKmersUsed);
                    }
                }

                log.info("");

                stretchNum++;
            }
        }

        for (String vid : viSeen.keySet()) {
            if (!viSeen.get(vid)) {
                evalTables.getTable("discoveryStats").increment("dummy", "fn");

                evalTables.getTable("variantStats").set(vid, "knownVariantId", vid);
                evalTables.getTable("variantStats").set(vid, "knownVariantEvent", vids.get(vid).denovo);
            }
        }

        evalTables.write(out);
    }
}
