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

    @Argument(fullName = "maxJunctionsAllowed", shortName = "j", doc = "Maximum number of junctions we'll permit ourselves to traverse")
    public Integer MAX_JUNCTIONS_ALLOWED = 10;

    @Argument(fullName="skipToKmer", shortName="s", doc="Skip processing to given kmer", required=false)
    public String KMER;

    @Output(fullName = "gout", shortName = "go", doc = "Graph out")
    public File gout;

    @Output
    public PrintStream out;

    @Output(fullName = "evalOut", shortName = "eout", doc = "Eval out")
    public PrintStream eout;

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

            log.info("  {}", gout.getAbsolutePath() + "." + prefix + ".pdf");

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
                    attrs.put("shape", withText ? "box" : "circle");
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
                        attrs.put("weight", g.getEdgeWeight(e));
                        attrs.put("label", g.getEdgeWeight(e));

                        o.println(indent + "\"" + s + "\" -> \"" + t + "\" [ " + formatAttributes(attrs) + " ];");
                    }
                }
            }

            o.println("}");

            o.close();

            if (withPdf) {
                Runtime.getRuntime().exec("dot -Tpdf -o" + p.getAbsolutePath() + " " + f.getAbsolutePath());
            }
        } catch (FileNotFoundException e) {
            throw new IndianaException("File not found", e);
        } catch (IOException e) {
            throw new IndianaException("IO exception", e);
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

    private DirectedGraph<AnnotatedVertex, AnnotatedEdge> simplifyGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, boolean combine) {
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

                while (nextVertex != null && b.outDegreeOf(thisVertex) == 1 && b.inDegreeOf(nextVertex) == 1 && (combine || thisVertex.isNovel() == nextVertex.isNovel())) {
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
                String knownAlt = bestVi.alt != null ? bestVi.alt : "";
                String ref = gvc.getAttributeAsString(color, "parentalAllele");
                String alt = gvc.getAttributeAsString(color, "childAllele");

                gvc.attribute(color, "isKnownVariant", true);
                gvc.attribute(color, "variantId", bestVi.variantId);
                gvc.attribute(color, "allelesMatch", (knownRef.equals(ref) && knownAlt.equals(alt)) || bestVi.denovo.equals(gvc.getAttributeAsString(color, "event")));
                gvc.attribute(color, "eventsMatch", bestVi.denovo.equals(gvc.getAttributeAsString(color, "event")));
                gvc.attribute(color, "knownRef", knownRef);
                gvc.attribute(color, "knownAlt", knownAlt);
            } else {
                gvc.attribute(color, "isKnownVariant", false);
                gvc.attribute(color, "variantId", "unknown");
                gvc.attribute(color, "allelesMatch", false);
                gvc.attribute(color, "eventsMatch", false);
                gvc.attribute(color, "knownRef", "unknown");
                gvc.attribute(color, "knownAlt", "unknown");
            }
        }
    }

    /*
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
    */

    private String recordToString(String sk, CortexRecord cr) {
        String kmer = cr.getKmerAsString();
        String cov = "";
        String ed = "";

        boolean fw = sk.equals(kmer);

        if (!fw) {
            kmer = SequenceUtils.reverseComplement(kmer);
        }

        for (int coverage : cr.getCoverages()) {
            cov += " " + coverage;
        }

        for (String edge : cr.getEdgeAsStrings()) {
            ed += " " + (fw ? edge : SequenceUtils.reverseComplement(edge));
        }

        return kmer + " " + cov + " " + ed;
    }

    @Override
    public void execute() {
        log.info("Loading reference indices for fast kmer lookup...");
        KmerLookup kl1 = new KmerLookup(REF1);
        KmerLookup kl2 = new KmerLookup(REF2);

        Map<CortexKmer, VariantInfo> vis = GenotypeGraphUtils.loadNovelKmerMap(NOVEL_KMER_MAP, BED);
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

        if (KMER != null) {
            novelKmers.clear();
            novelKmers.put(new CortexKmer(KMER), true);
        }

        int variantsMissingKmers = 0;

        log.info("Genotyping novel kmer stretches in graph...");
        Set<GraphicalVariantContext> gvcs = new LinkedHashSet<GraphicalVariantContext>();
        for (CortexKmer novelKmer : novelKmers.keySet()) {
            if (novelKmers.get(novelKmer)) {
                // Walk the graph left and right of novelKmer and extract a novel stretch
                String stretch;
                if (GRAPH_RAW == null) {
                    stretch = CortexUtils.getSeededStretch(GRAPH, novelKmer.getKmerAsString(), 0, AGGRESSIVE);
                } else {
                    stretch = CortexUtils.getSeededStretch(GRAPH, GRAPH_RAW, novelKmer.getKmerAsString(), 0, AGGRESSIVE);
                }

                log.info("  stretch {}: {} bp", stretchNum, stretch.length());
                log.info("    novel kmer: {}", novelKmer);
                log.info("    sequence:");
                log.info("    - 0: {}", SequenceUtils.truncate(stretch, 100));

                // Construct GVC
                GraphicalVariantContext gvc = new GraphicalVariantContext()
                        .attribute(0, "stretch", stretch)
                        .attribute(0, "stretchNum", stretchNum)
                        .attribute(0, "stretchLength", stretch.length())
                        .attribute(0, "novelKmersTotal", novelKmers.size());

                // Fetch the local subgraph context from disk
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> ag = GenotypeGraphUtils.dfsGraph(stretch, GRAPH, GRAPH_RAW, null, AGGRESSIVE, novelKmers);
                log.info("    subgraph : {} vertices, {} edges", ag.vertexSet().size(), ag.edgeSet().size());

                //log.debug("Graph printed");
                //printGraph(simplifyGraph(ag, false), "call" + String.format("%04d", stretchNum), false, true);

                if (true) {
                    // Extract parental stretches
                    PathInfo p1 = GenotypeGraphUtils.computeBestMinWeightPath(GRAPH, GRAPH_RAW, ag, 1, stretch, novelKmers);
                    PathInfo p2 = GenotypeGraphUtils.computeBestMinWeightPath(GRAPH, GRAPH_RAW, ag, 2, stretch, novelKmers);

                    log.info("    paths:");
                    log.info("    - 1: {} ({} bp)", SequenceUtils.truncate(p1.parent, 100), p1.parent.length());
                    log.info("      c: {} ({} bp)", SequenceUtils.truncate(p1.child, 100), p1.child.length());
                    log.info("    - 2: {} ({} bp)", SequenceUtils.truncate(p2.parent, 100), p2.parent.length());
                    log.info("      c: {} ({} bp)", SequenceUtils.truncate(p2.child, 100), p2.child.length());

                    // Call variants
                    gvc.add(GenotypeGraphUtils.callVariant(GRAPH, GRAPH_RAW, p1, 1, stretch, novelKmers, kl1));
                    gvc.add(GenotypeGraphUtils.callVariant(GRAPH, GRAPH_RAW, p2, 2, stretch, novelKmers, kl2));

                    log.info("    variants:");
                    log.info("    - 1: {} {} ({} bp)", gvc.getAttributeAsString(1, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(1, "parentalAllele"), 70), gvc.getAttributeAsString(1, "parentalAllele").length());
                    log.info("      c: {} {} ({} bp)", gvc.getAttributeAsString(1, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(1, "childAllele"), 70), gvc.getAttributeAsString(1, "childAllele").length());
                    log.info("    - 2: {} {} ({} bp)", gvc.getAttributeAsString(2, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(2, "parentalAllele"), 70), gvc.getAttributeAsString(2, "parentalAllele").length());
                    log.info("      c: {} {} ({} bp)", gvc.getAttributeAsString(2, "event"), SequenceUtils.truncate(gvc.getAttributeAsString(2, "childAllele"), 70), gvc.getAttributeAsString(2, "childAllele").length());

                    // Finalize into a single call
                    GenotypeGraphUtils.chooseVariant(gvc);

                    // Show alignment
                    log.info("    alignment:");
                    log.info("    - novel stretch: {}", gvc.getAttribute(0, "novelStretchAlignment"));
                    log.info("    - parental path: {}", gvc.getAttribute(0, "parentalPathAlignment"));
                }

                // See how many novel kmers we've used up
                int novelKmersUsed = 0;

                for (AnnotatedVertex av : ag.vertexSet()) {
                    CortexKmer ck = new CortexKmer(av.getKmer());

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

                        log.info("    - {}: background {}, isKnownVariant {}, allelesMatch {}, eventsMatch {}, variantId {} {} {} {}",
                                c,
                                gvc.getAttributeAsInt(c, "haplotypeBackground"),
                                gvc.getAttribute(c, "isKnownVariant"),
                                gvc.getAttribute(c, "allelesMatch"),
                                gvc.getAttribute(c, "eventsMatch"),
                                gvc.getAttributeAsString(c, "variantId"),
                                event,
                                gvc.getAttributeAsString(c, "knownRef"),
                                gvc.getAttributeAsString(c, "knownAlt")
                        );

                        /*
                        if (vi != null) {
                            String refSeq = vi.leftFlank + vi.ref + vi.rightFlank;
                            String altSeq = vi.leftFlank + (vi.alt == null ? "" : vi.alt) + vi.rightFlank;

                            boolean isFwd = true;
                            boolean isMissingKmers = false;

                            for (int i = 0; i <= refSeq.length() - GRAPH.getKmerSize(); i++) {
                                String fw = refSeq.substring(i, i + GRAPH.getKmerSize());
                                String rc = SequenceUtils.reverseComplement(fw);
                                CortexRecord cr = GRAPH.findRecord(new CortexKmer(fw));

                                if (ag.containsVertex(new AnnotatedVertex(fw)) || ag.containsVertex(new AnnotatedVertex(fw, true))) {
                                    log.info("    - ref {}/{}: fw {} {}", i, refSeq.length() - GRAPH.getKmerSize(), fw, recordToString(fw, cr));

                                    isFwd = true;
                                } else if (ag.containsVertex(new AnnotatedVertex(rc)) || ag.containsVertex(new AnnotatedVertex(rc, true))) {
                                    log.info("    - ref {}/{}: rc {} {}", i, refSeq.length() - GRAPH.getKmerSize(), rc, recordToString(rc, cr));
                                    isFwd = false;
                                } else {
                                    log.info("    - ref {}/{}: ?? {} {}", i, refSeq.length() - GRAPH.getKmerSize(), isFwd ? fw : rc, recordToString(isFwd ? fw : rc, cr));
                                    isMissingKmers = true;
                                }
                            }


                            for (int i = 0; i <= altSeq.length() - GRAPH.getKmerSize(); i++) {
                                String fw = altSeq.substring(i, i + GRAPH.getKmerSize());
                                String rc = SequenceUtils.reverseComplement(fw);

                                //AnnotatedVertex afw = new AnnotatedVertex(fw);
                                //AnnotatedVertex arc = new AnnotatedVertex(SequenceUtils.reverseComplement(fw));

                                if (ag.containsVertex(new AnnotatedVertex(fw)) || ag.containsVertex(new AnnotatedVertex(fw, true))) {
                                    //log.info("    - {}: fw {}", i, fw);
                                    isFwd = true;
                                } else if (ag.containsVertex(new AnnotatedVertex(rc)) || ag.containsVertex(new AnnotatedVertex(rc, true))) {
                                    //log.info("    - {}: rc {}", i, rc);
                                    isFwd = false;
                                } else {
                                    log.info("    - alt {}/{}: ?? {}", i, altSeq.length() - GRAPH.getKmerSize(), isFwd ? fw : rc);
                                    isMissingKmers = true;
                                }
                            }

                            if (isMissingKmers) {
                                variantsMissingKmers++;
                            }
                        }

                        break;
                        */
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
                        evalTables.getTable("variantStats").set(pk, "seedKmer", novelKmer.getKmerAsString());

                        viSeen.put(vi.variantId, true);

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
                        evalTables.getTable("variantStats").set(pk, "seedKmer", novelKmer.getKmerAsString());
                    }
                }

                log.info("");

                gvcs.add(gvc);
                stretchNum++;
            }
        }

        log.info("Variants missing kmers: {}", variantsMissingKmers);

        for (String vid : viSeen.keySet()) {
            if (!viSeen.get(vid)) {
                evalTables.getTable("discoveryStats").increment("dummy", "fn");

                evalTables.getTable("variantStats").set(vid, "knownVariantId", vid);
                evalTables.getTable("variantStats").set(vid, "knownVariantEvent", vids.get(vid).denovo);

                String pk = vid + ".none";

                VariantInfo vi = vids.get(vid);

                int refLength = vi.ref != null ? vi.ref.length() : 0;
                int altLength = vi.alt != null ? vi.alt.length() : 0;

                evalTables.getTable("variantStats").set(pk, "knownVariantId", vid);
                evalTables.getTable("variantStats").set(pk, "knownVariantEvent", vids.get(vid).denovo);
                evalTables.getTable("variantStats").set(pk, "knownVariantLength", Math.abs(refLength - altLength));
                evalTables.getTable("variantStats").set(pk, "variantId", "none");
                evalTables.getTable("variantStats").set(pk, "variantEvent", "none");
                evalTables.getTable("variantStats").set(pk, "variantLength", 0);
                evalTables.getTable("variantStats").set(pk, "novelKmersUsed", 0);
                evalTables.getTable("variantStats").set(pk, "seedKmer", "none");
            }
        }

        if (BED != null) {
            evalTables.write(out);
        }
    }
}
