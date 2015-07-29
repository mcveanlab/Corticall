package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.DirectedGraph;
import org.jgrapht.GraphPath;
import org.jgrapht.Graphs;
import org.jgrapht.alg.DijkstraShortestPath;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.EdgeReversedGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.exact.ExactLookup;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public class GenotypeGraph extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="graphRaw", shortName="r", doc="Graph (raw)", required=false)
    public CortexGraph GRAPH_RAW;

    @Argument(fullName="novelGraph", shortName="n", doc="Graph of novel kmers")
    public CortexGraph NOVEL;

    @Argument(fullName="ref1", shortName="r1", doc="Fasta file for first parent")
    public File REF1;

    @Argument(fullName="ref2", shortName="r2", doc="Fasta file for second parent")
    public File REF2;

    @Argument(fullName="bed", shortName="b", doc="Bed file describing variants", required=false)
    public File BED;

    @Argument(fullName="novelKmerMap", shortName="m", doc="Novel kmer map", required=false)
    public File NOVEL_KMER_MAP;

    @Output(fullName="gout", shortName="go", doc="Graph out")
    public File gout;

    @Output
    public PrintStream out;

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

                    if (kv[0].equals("id")) { vi.variantId = kv[1]; }
                    if (kv[0].equals("type")) { vi.type = kv[1]; }
                    if (kv[0].equals("denovo")) { vi.denovo = kv[1]; }
                    if (kv[0].equals("nahr")) { vi.nahr = kv[1]; }
                    if (kv[0].equals("gcindex")) { vi.gcindex = Integer.valueOf(kv[1]); }
                    if (kv[0].equals("ref")) { vi.ref = kv[1]; }
                    if (kv[0].equals("alt")) { vi.alt = kv[1]; }
                    if (kv[0].equals("left")) { vi.leftFlank = kv[1]; }
                    if (kv[0].equals("right")) { vi.rightFlank = kv[1]; }
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

            String[] colors = new String[] { "red", "blue", "green" };

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

    private void addGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, DirectedGraph<String, DefaultEdge> g, int color, Map<CortexKmer, Boolean> novelKmers) {
        for (String v : g.vertexSet()) {
            AnnotatedVertex av = new AnnotatedVertex(v, novelKmers.containsKey(new CortexKmer(v)));

            a.addVertex(av);
        }

        for (DefaultEdge e : g.edgeSet()) {
            String s0 = g.getEdgeSource(e);
            String s1 = g.getEdgeTarget(e);

            CortexKmer ck0 = new CortexKmer(s0);
            CortexKmer ck1 = new CortexKmer(s1);

            AnnotatedVertex a0 = new AnnotatedVertex(s0, novelKmers.containsKey(ck0));
            AnnotatedVertex a1 = new AnnotatedVertex(s1, novelKmers.containsKey(ck1));

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

    private Pair<String, String> computeBestMinWeightPath(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, int color, String stretch, Map<CortexKmer, Boolean> novelKmers) {
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> b = copyGraph(a);

        AnnotateStartsAndEnds annotateStartsAndEnds = new AnnotateStartsAndEnds(color, stretch, novelKmers, b).invoke();
        Set<AnnotatedVertex> candidateStarts = annotateStartsAndEnds.getCandidateStarts();
        Set<AnnotatedVertex> candidateEnds = annotateStartsAndEnds.getCandidateEnds();

        DirectedGraph<AnnotatedVertex, AnnotatedEdge> b0 = removeOtherColors(b, 0);
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> bc = removeOtherColors(b, color);

        List<String> novelStretches = getNovelStretch(stretch, novelKmers);

        double minPl0 = Double.MAX_VALUE, minPlc = Double.MAX_VALUE;
        String minLp0 = "", minLpc = "";

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
                        }
                    }
                }
            }
        }

        return new Pair<String, String>(minLp0, minLpc);
    }

    private enum KmerOrigin { MOTHER, FATHER, BOTH, NONE };

    private void align(String stretch, KmerLookup kl) {
        KmerOrigin[] ko = new KmerOrigin[stretch.length() - GRAPH.getKmerSize() + 1];

        List<String> sks = new ArrayList<String>();
        for (int i = 0; i <= stretch.length() - GRAPH.getKmerSize(); i++) {
            String sk = stretch.substring(i, i + GRAPH.getKmerSize());
            CortexKmer ck = new CortexKmer(sk);
            CortexRecord cr = GRAPH.findRecord(ck);

            sks.add(sk);

            //log.info("sk: {} {}", sk, kl.find(sk));

            //List<Set<Interval>> interval = kl.find(sk);

            if (cr != null) {
                int cov1 = cr.getCoverage(1);
                int cov2 = cr.getCoverage(2);

                if (cov1 > 0 && cov2 > 0) {
                    ko[i] = KmerOrigin.BOTH;
                } else if (cov1 >  0 && cov2 == 0) {
                    ko[i] = KmerOrigin.MOTHER;
                } else if (cov1 == 0 && cov2 >  0) {
                    ko[i] = KmerOrigin.FATHER;
                } else {
                    ko[i] = KmerOrigin.NONE;
                }
            } else {
                ko[i] = KmerOrigin.NONE;
            }
        }

        int start = 0, end = 0;
        KmerOrigin currentOrigin = null;
        int stretchNum = 0;
        Set<String> pieces = new LinkedHashSet<String>();

        for (int i = 0; i < ko.length; i++) {
            String sk = stretch.substring(i, i + GRAPH.getKmerSize());

            if ((ko[i] == KmerOrigin.MOTHER || ko[i] == KmerOrigin.FATHER) && currentOrigin != null && currentOrigin == KmerOrigin.BOTH) {
                currentOrigin = ko[i];
            }

            if (currentOrigin == null && ko[i] != KmerOrigin.NONE) {
                currentOrigin = ko[i];
                start = i;
                end = i;
            } else if (ko[i] == currentOrigin || ko[i] == KmerOrigin.BOTH) {
                if (currentOrigin == KmerOrigin.BOTH && ko[i] != currentOrigin) {
                    currentOrigin = ko[i];
                }

                end = i;
            } else if (ko[i] != currentOrigin && ko[i] != KmerOrigin.BOTH && currentOrigin != null) {
                StringBuilder sb = new StringBuilder(sks.get(start));

                for (int j = start + 1; j <= end; j++) {
                    String skj = sks.get(j);
                    sb.append(skj.charAt(skj.length() - 1));
                }

                pieces.add(sb.toString());

                stretchNum++;
                start = i;
                end = i;
                currentOrigin = (ko[i] == KmerOrigin.NONE) ? null : ko[i];
            }

            //log.info("{} {} {} {} {}-{}", i, sk, ko[i], stretchNum, start, end);
        }

        if (currentOrigin != null && currentOrigin != KmerOrigin.NONE) {
            StringBuilder sb = new StringBuilder(sks.get(start));

            for (int j = start + 1; j <= end; j++) {
                String skj = sks.get(j);
                sb.append(skj.charAt(skj.length() - 1));
            }

            pieces.add(sb.toString());
        }

        for (String piece : pieces) {
            log.info("piece: {} {}", kl.find(piece), piece);
        }
    }

    private VariantContext callVariant(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, int color, String stretch, Map<CortexKmer, Boolean> novelKmers, KmerLookup kl) {
        Pair<String, String> p = computeBestMinWeightPath(a, color, stretch, novelKmers);

        int kmerSize = GRAPH.getKmerSize();
        String startFw, endFw;
        if (p.getSecond().length() > GRAPH.getKmerSize()) {
            startFw = p.getSecond().substring(0, kmerSize);
            endFw = p.getSecond().substring(p.getSecond().length() - kmerSize, p.getSecond().length());
        } else {
            startFw = null;
            endFw = null;

            for (int i = 0; i <= stretch.length() - kmerSize; i++) {
                String sk = stretch.substring(i, i + kmerSize);
                CortexKmer ck = new CortexKmer(sk);
                CortexRecord cr = GRAPH.findRecord(ck);

                if (cr != null && cr.getCoverage(color) > 0) { startFw = sk; }
            }

            for (int i = stretch.length() - kmerSize; i >= 0; i--) {
                String sk = stretch.substring(i, i + kmerSize);
                CortexKmer ck = new CortexKmer(sk);
                CortexRecord cr = GRAPH.findRecord(ck);

                if (cr != null && cr.getCoverage(color) > 0) { endFw = sk; }
            }
        }

        if (startFw != null && endFw != null) {
            String contextLeft = CortexUtils.getSeededStretchLeft(GRAPH, startFw, color);
            String contextRight = CortexUtils.getSeededStretchRight(GRAPH, endFw, color);
            String contextFull = contextLeft + p.getSecond() + contextRight;

            List<Set<Interval>> intervalLeft = kl.find(contextLeft);
            List<Set<Interval>> intervalRight = kl.find(contextRight);
            List<Set<Interval>> intervalFull = kl.find(contextFull);

            log.info("    context:");
            log.info("    - {} {} {} {}", color, intervalLeft, contextLeft.length(), contextLeft);
            log.info("    - {} {} {} {}", color, intervalRight, contextRight.length(), contextRight);
            log.info("    - {} {} {} {}", color, intervalFull, contextFull.length(), contextFull);

            if (contextFull.length() > GRAPH.getKmerSize()) {
                align(contextFull, kl);
            } else {
                align(contextLeft + stretch + contextRight, kl);
            }
        }

        int s, e0 = p.getFirst().length() - 1, e1 = p.getSecond().length() - 1;

        for (s = 0; s < (p.getFirst().length() < p.getSecond().length() ? p.getFirst().length() : p.getSecond().length()) && p.getFirst().charAt(s) == p.getSecond().charAt(s); s++) {}

        while (e0 > s && e1 > s && p.getFirst().charAt(e0) == p.getSecond().charAt(e1)) {
            e0--;
            e1--;
        }

        String parentalAllele = p.getSecond() == null || p.getSecond().equals("") || p.getFirst().equals(p.getSecond()) ? "A" : p.getSecond().substring(s, e1 + 1);
        String childAllele = p.getFirst() == null || p.getFirst().equals("") || p.getFirst().equals(p.getSecond()) ? "N" : p.getFirst().substring(s, e0 + 1);

        int e = s + parentalAllele.length() - 1;

        VariantContextBuilder vcb = new VariantContextBuilder()
                .chr("unknown")
                .start(s)
                .stop(e)
                .alleles(parentalAllele, childAllele)
                .attribute("NOVEL_STRETCH", stretch)
                .attribute("CHILD_STRETCH", p.getFirst())
                .attribute("PARENT_STRETCH", p.getSecond())
                .attribute("CHILD_COLOR", 0)
                .attribute("PARENT_COLOR", color)
                .source(String.valueOf(color));

        if (childAllele.equals("N")) {
            vcb.filter("TRAVERSAL_INCOMPLETE");
        }

        return vcb.make();
    }

    private DirectedGraph<AnnotatedVertex, AnnotatedEdge> loadLocalGraph(Map<CortexKmer, Boolean> novelKmers, CortexKmer novelKmer, String stretch) {
        // first, explore each color and bring the local subgraphs into memory
        DirectedGraph<String, DefaultEdge> sg0 = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
        DirectedGraph<String, DefaultEdge> sg1 = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
        DirectedGraph<String, DefaultEdge> sg2 = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);

        for (int i = 0; i <= stretch.length() - novelKmer.length(); i++) {
            String kmer = stretch.substring(i, i + novelKmer.length());

            if (!sg0.containsVertex(kmer)) {
                Graphs.addGraph(sg0, CortexUtils.dfs(GRAPH, kmer, 0, null, new AbstractTraversalStopper() {
                    @Override
                    public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<String, DefaultEdge> g, int depth) {
                        return cr.getCoverage(1) > 0 || cr.getCoverage(2) > 0;
                    }

                    @Override
                    public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<String, DefaultEdge> g, int junctions) {
                        //return junctions >= maxJunctionsAllowed();
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
            for (String kmer : sg0.vertexSet()) {
                DirectedGraph<String, DefaultEdge> sg = (c == 1) ? sg1 : sg2;

                if (!sg.containsVertex(kmer)) {
                    TraversalStopper stopper = new AbstractTraversalStopper() {
                        @Override
                        public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<String, DefaultEdge> g, int junctions) {
                            String fw = cr.getKmerAsString();
                            String rc = SequenceUtils.reverseComplement(fw);

                            if (g.containsVertex(fw) || g.containsVertex(rc)) {
                                if (junctions < distanceToGoal) {
                                    distanceToGoal = junctions;
                                }

                                return true;
                            }

                            return false;
                        }

                        @Override
                        public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<String, DefaultEdge> g, int junctions) {
                            return junctions >= maxJunctionsAllowed();
                        }

                        @Override
                        public int maxJunctionsAllowed() {
                            return 5;
                        }
                    };

                    Graphs.addGraph(sg, CortexUtils.dfs(GRAPH, kmer, c, sg0, stopper));
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

            String fk = stretch.substring(1, kmerLength + 1);
            AnnotatedVertex afk = new AnnotatedVertex(fk, novelKmers.containsKey(new CortexKmer(fk)));

            String lk = stretch.substring(stretch.length() - kmerLength - 1, stretch.length() - 1);
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

    private boolean evalVariant(VariantContext vc, Map<CortexKmer, VariantInfo> vis, String stretch) {
        Set<VariantInfo> relevantVis = new HashSet<VariantInfo>();
        int kmerSize = vis.keySet().iterator().next().length();

        for (int i = 0; i <= stretch.length() - kmerSize; i++) {
            CortexKmer ck = new CortexKmer(stretch.substring(i, i + kmerSize));

            if (vis.containsKey(ck)) {
                relevantVis.add(vis.get(ck));
            }
        }

        if (relevantVis.size() > 0) {
            VariantInfo vi = relevantVis.iterator().next();

            String ref = vc.getReference().getBaseString();
            String alt = vc.getAlternateAllele(0).getBaseString();

            String refStretch = vc.getAttributeAsString("PARENT_STRETCH", ref);
            String altStretch = vc.getAttributeAsString("CHILD_STRETCH", alt);

            int pos = vc.getStart();
            int refLength = vc.getReference().length();
            int altLength = vc.getAlternateAllele(0).length();
            boolean found = false;

            while (pos >= 0 && pos + refLength < refStretch.length() && pos + altLength < altStretch.length()) {
                String refFw = refStretch.substring(pos, pos + refLength);
                String refRc = SequenceUtils.reverseComplement(refFw);

                String altFw = altStretch.substring(pos, pos + altLength);
                String altRc = SequenceUtils.reverseComplement(altFw);

                if (vi.ref.equals(refFw) && vi.alt.equals(altFw)) {
                    ref = refFw;
                    alt = altFw;
                    found = true;
                    break;
                } else if (vi.ref.equals(refRc) && vi.alt.equals(altRc)) {
                    ref = refRc;
                    alt = altRc;
                    found = true;
                    break;
                }

                pos--;
            }

            if (!found) {
                pos = vc.getStart();

                while (pos >= 0 && pos + refLength < refStretch.length() && pos + altLength < altStretch.length()) {
                    String refFw = refStretch.substring(pos, pos + refLength);
                    String refRc = SequenceUtils.reverseComplement(refFw);

                    String altFw = altStretch.substring(pos, pos + altLength);
                    String altRc = SequenceUtils.reverseComplement(altFw);

                    if (vi.ref.equals(refFw) && vi.alt.equals(altFw)) {
                        ref = refFw;
                        alt = altFw;
                        found = true;
                        break;
                    } else if (vi.ref.equals(refRc) && vi.alt.equals(altRc)) {
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

            log.info("    - vi: {}", vi);
            log.info("        known ref: {}", knownRef);
            log.info("        known alt: {}", knownAlt);
            log.info("        called ref: {}", ref);
            log.info("        called alt: {}", alt);

            log.info("    - matches: {} ({} {} {} {})", knownRef.equals(ref) && knownAlt.equals(alt), ref, alt, knownRef, knownAlt);

            //return (knownRef.equals(ref) && knownAlt.equals(alt)) || vi.denovo.equals("GC") || vi.denovo.equals("NAHR");
            return (knownRef.equals(ref) && knownAlt.equals(alt)) || vi.denovo.equals("GC");
        }

        log.info("    - matches: {}", false);

        return false;
    }

    @Override
    public void execute() {
        log.info("Loading reference indices for fast kmer lookup...");
        KmerLookup kl1 = new KmerLookup(REF1);
        KmerLookup kl2 = new KmerLookup(REF2);

        Map<CortexKmer, VariantInfo> vis = loadNovelKmerMap();

        Map<CortexKmer, Boolean> novelKmers = new HashMap<CortexKmer, Boolean>();
        for (CortexRecord cr : NOVEL) {
            novelKmers.put(new CortexKmer(cr.getKmerAsString()), true);
        }

        //TableWriter tw = new TableWriter(out);
        int totalNovelKmersUsed = 0;
        int stretchNum = 0;
        int numMatches = 0;

        log.info("Genotyping novel kmer stretches in graph...");
        for (CortexKmer novelKmer : novelKmers.keySet()) {
            if (novelKmers.get(novelKmer)) {
                int novelKmersUsed = 0;

                // Walk the graph left and right of novelKmer and extract a novel stretch
                String stretch = CortexUtils.getSeededStretch(GRAPH, novelKmer.getKmerAsString(), 0);
                log.info("  stretch : {} ({} bp)", stretchNum, stretch.length());

                // Fetch the local subgraph context from disk
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> ag = loadLocalGraph(novelKmers, novelKmer, stretch);
                log.info("    subgraph: {} vertices, {} edges", ag.vertexSet().size(), ag.edgeSet().size());

                // See how many novel kmers we've used up
                for (AnnotatedVertex ak : ag.vertexSet()) {
                    CortexKmer ck = new CortexKmer(ak.getKmer());

                    if (novelKmers.containsKey(ck) && novelKmers.get(ck)) {
                        totalNovelKmersUsed++;
                        novelKmersUsed++;
                        novelKmers.put(ck, false);
                    }
                }

                log.info("    novelty : novelKmersUsed: {}, totalNovelKmersUsed: {}/{}", novelKmersUsed, totalNovelKmersUsed, novelKmers.size());

                // Extract stretches
                Pair<String, String> p1 = computeBestMinWeightPath(ag, 1, stretch, novelKmers);
                Pair<String, String> p2 = computeBestMinWeightPath(ag, 2, stretch, novelKmers);

                log.info("    stretches:");
                log.info("    - s1: {}", p1.getSecond());
                log.info("          {}", p1.getFirst());
                log.info("    - s2: {}", p2.getSecond());
                log.info("          {}", p2.getFirst());

                // Call variants
                VariantContext vc1 = callVariant(ag, 1, stretch, novelKmers, kl1);
                VariantContext vc2 = callVariant(ag, 2, stretch, novelKmers, kl2);

                log.info("    variants:");
                log.info("    - vc1: {}", vc1);
                log.info("    - vc2: {}", vc2);

                // Evaluate variants
                if (BED != null) {
                    log.info("    evaluate:");
                    boolean m1 = evalVariant(vc1, vis, stretch);
                    boolean m2 = evalVariant(vc2, vis, stretch);

                    printGraph(simplifyGraph(ag), "debug", false, false);
                    printGraph(ag, "debugFull", false, false);

                    if (m1 || m2) {
                        log.info("    - match!");
                        numMatches++;
                    } else {
                        log.info("    - no match :(");

                        for (AnnotatedVertex v : ag.vertexSet()) {
                            if (v.isNovel() && (ag.inDegreeOf(v) == 0 || ag.outDegreeOf(v) == 0)) {
                                log.info("weird: {}", v);
                            }
                        }

                        boolean b1 = evalVariant(vc1, vis, stretch);
                        boolean b2 = evalVariant(vc2, vis, stretch);

                        Pair<String, String> d1 = computeBestMinWeightPath(ag, 1, stretch, novelKmers);
                        Pair<String, String> d2 = computeBestMinWeightPath(ag, 2, stretch, novelKmers);
                    }
                }

                stretchNum++;
            }
        }

        log.info("Num stretches: {}", stretchNum);
        log.info("Num matches: {}", numMatches);
    }
}
