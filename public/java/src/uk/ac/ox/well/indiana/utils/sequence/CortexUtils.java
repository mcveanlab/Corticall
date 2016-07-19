package uk.ac.ox.well.indiana.utils.sequence;

import com.google.common.base.Joiner;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultDirectedGraph;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalStopper;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;

/**
 * A set of utilities for dealing with Cortex graphs and records
 */
public class CortexUtils {
    /**
     * Private constructor - this class cannot be instantiated!
     */
    private CortexUtils() {}

    public static String getNextKmer(CortexGraph cg, String kmer, int color, boolean beAggressive) {
        CortexKmer ck = new CortexKmer(kmer);
        CortexRecord cr = cg.findRecord(ck);

        if (cr != null) {
            Collection<String> outEdges = cr.getOutEdgesAsStrings(color);

            if (ck.isFlipped()) {
                outEdges = cr.getInEdgesComplementAsStrings(color);
            }

            if (outEdges.size() == 1) {
                String outEdge = outEdges.iterator().next();

                return kmer.substring(1, kmer.length()) + outEdge;
            } else if (beAggressive) {
                Set<String> inferredKmers = new HashSet<>();
                Set<String> novelKmers = new HashSet<>();

                for (String outEdge : Arrays.asList("A", "C", "G", "T")) {
                    String outKmer = kmer.substring(1, kmer.length()) + outEdge;

                    if (hasKmer(cg, outKmer, color)) {
                        inferredKmers.add(outKmer);

                        if (isNovelKmer(cg, outKmer, color)) {
                            novelKmers.add(outKmer);
                        }
                    }
                }

                if (novelKmers.size() == 1) {
                    return novelKmers.iterator().next();
                } else if (inferredKmers.size() == 1) {
                    return inferredKmers.iterator().next();
                }
            }
        }

        return null;
    }

    public static String getNextKmer(CortexGraph clean, CortexGraph dirty, String kmer, int color, boolean beAggressive) {
        String nextClean = getNextKmer(clean, kmer, color, beAggressive);
        String nextDirty = null;

        if (dirty != null) {
            nextDirty = getNextKmer(dirty, kmer, color, beAggressive);
        }

        if (nextClean != null) { return nextClean; }
        else if (nextDirty != null) { return nextDirty; }

        return null;
    }

    public static Set<String> getNextKmers(CortexGraph cg, String kmer, int color) {
        CortexKmer ck = new CortexKmer(kmer);
        CortexRecord cr = cg.findRecord(ck);
        Set<String> nextKmers = new HashSet<>();

        if (cr != null) {
            Collection<String> outEdges = cr.getOutEdgesAsStrings(color);

            if (ck.isFlipped()) {
                outEdges = cr.getInEdgesComplementAsStrings(color);
            }

            for (String outEdge : outEdges) {
                nextKmers.add(kmer.substring(1, kmer.length()) + outEdge);
            }
        }

        return nextKmers;
    }

    public static Set<String> getNextKmers(CortexGraph clean, CortexGraph dirty, String kmer, int color) {
        Set<String> nextKmers = getNextKmers(clean, kmer, color);

        if (dirty != null && (nextKmers == null || nextKmers.isEmpty())) {
            nextKmers = getNextKmers(dirty, kmer, color);
        }

        return nextKmers;
    }

    private static Pair<CortexGraph, CortexRecord> getGraphRecordPair(CortexGraph clean, CortexGraph dirty, String kmer) {
        CortexKmer ck = new CortexKmer(kmer);

        CortexRecord crclean = clean == null ? null : clean.findRecord(ck);
        CortexRecord crdirty = dirty == null ? null : dirty.findRecord(ck);

        CortexRecord cr = crclean == null && crdirty != null ? crdirty : crclean;
        CortexGraph cg = clean == null && dirty != null ? dirty : clean;

        return new Pair<>(cg, cr);
    }

    public static Map<Integer, Set<String>> getNextKmers(CortexGraph clean, CortexGraph dirty, String kmer) {
        Pair<CortexGraph, CortexRecord> cgcr = getGraphRecordPair(clean, dirty, kmer);

        return getNextKmers(cgcr, kmer);
    }

    public static Map<Integer, Set<String>> getNextKmers(Pair<CortexGraph, CortexRecord> cgcr, String kmer) {
        CortexKmer ck = new CortexKmer(kmer);

        CortexGraph cg = cgcr.getFirst();
        CortexRecord cr = cgcr.getSecond();

        Map<Integer, Set<String>> nextKmersAllColors = new HashMap<>();

        if (cg != null) {
            for (int c = 0; c < cg.getNumColors(); c++) {
                Set<String> nextKmers = new HashSet<>();

                if (cr != null) {
                    Collection<String> outEdges = (ck.isFlipped()) ? cr.getInEdgesComplementAsStrings(c) : cr.getOutEdgesAsStrings(c);

                    for (String outEdge : outEdges) {
                        nextKmers.add(kmer.substring(1, kmer.length()) + outEdge);
                    }
                }

                nextKmersAllColors.put(c, nextKmers);
            }
        }

        return nextKmersAllColors;
    }

    public static String getPrevKmer(CortexGraph cg, String kmer, int color, boolean beAggressive) {
        CortexKmer ck = new CortexKmer(kmer);
        CortexRecord cr = cg.findRecord(ck);

        if (cr != null) {
            Collection<String> inEdges = cr.getInEdgesAsStrings(color);

            if (ck.isFlipped()) {
                inEdges = cr.getOutEdgesComplementAsStrings(color);
            }

            if (inEdges.size() == 1) {
                String inEdge = inEdges.iterator().next();

                return inEdge + kmer.substring(0, kmer.length() - 1);
            } else if (beAggressive) {
                Set<String> inferredKmers = new HashSet<>();
                Set<String> novelKmers = new HashSet<>();

                for (String inEdge : Arrays.asList("A", "C", "G", "T")) {
                    String inKmer = inEdge + kmer.substring(0, kmer.length() - 1);

                    if (hasKmer(cg, inKmer, color)) {
                        inferredKmers.add(inKmer);

                        if (isNovelKmer(cg, inKmer, color)) {
                            novelKmers.add(inKmer);
                        }
                    }
                }

                if (novelKmers.size() == 1) {
                    return novelKmers.iterator().next();
                } else if (inferredKmers.size() == 1) {
                    return inferredKmers.iterator().next();
                }
            }
        }

        return null;
    }

    public static String getPrevKmer(CortexGraph clean, CortexGraph dirty, String kmer, int color, boolean beAggressive) {
        String prevClean = getPrevKmer(clean, kmer, color, beAggressive);
        String prevDirty = null;

        if (dirty != null) {
            prevDirty = getPrevKmer(dirty, kmer, color, beAggressive);
        }

        if (prevClean != null) { return prevClean; }
        else if (prevDirty != null) { return prevDirty; }

        return null;
    }

    public static Set<String> getPrevKmers(CortexGraph cg, String kmer, int color) {
        CortexKmer ck = new CortexKmer(kmer);
        CortexRecord cr = cg.findRecord(ck);
        Set<String> prevKmers = new HashSet<>();

        if (cr != null) {
            Collection<String> inEdges = cr.getInEdgesAsStrings(color);

            if (ck.isFlipped()) {
                inEdges = cr.getOutEdgesComplementAsStrings(color);
            }

            for (String inEdge : inEdges) {
                prevKmers.add(inEdge + kmer.substring(0, kmer.length() - 1));
            }
        }

        return prevKmers;
    }

    public static Set<String> getPrevKmers(CortexGraph clean, CortexGraph dirty, String kmer, int color) {
        Set<String> prevKmers = getPrevKmers(clean, kmer, color);

        if (dirty != null && (prevKmers == null || prevKmers.isEmpty())) {
            prevKmers = getPrevKmers(dirty, kmer, color);
        }

        return prevKmers;
    }

    public static Map<Integer, Set<String>> getPrevKmers(CortexGraph clean, CortexGraph dirty, String kmer) {
        Pair<CortexGraph, CortexRecord> cgcr = getGraphRecordPair(clean, dirty, kmer);

        return getPrevKmers(cgcr, kmer);
    }

    public static Map<Integer, Set<String>> getPrevKmers(Pair<CortexGraph, CortexRecord> cgcr, String kmer) {
        CortexKmer ck = new CortexKmer(kmer);

        CortexGraph cg = cgcr.getFirst();
        CortexRecord cr = cgcr.getSecond();

        Map<Integer, Set<String>> prevKmersAllColors = new HashMap<>();

        if (cg != null) {
            for (int c = 0; c < cg.getNumColors(); c++) {
                Set<String> prevKmers = new HashSet<>();

                if (cr != null) {
                    Collection<String> inEdges = (ck.isFlipped()) ? cr.getOutEdgesComplementAsStrings(c) : cr.getInEdgesAsStrings(c);

                    for (String inEdge : inEdges) {
                        prevKmers.add(inEdge + kmer.substring(0, kmer.length() - 1));
                    }
                }

                prevKmersAllColors.put(c, prevKmers);
            }
        }

        return prevKmersAllColors;
    }

    public static String getSeededStretchLeft(CortexGraph cg, String kmer, int color, boolean beAggressive) {
        String tk = kmer;
        StringBuilder stretchBuilder = new StringBuilder();

        Set<String> usedKmers = new HashSet<>();

        String pk;
        while ((pk = CortexUtils.getPrevKmer(cg, tk, color, beAggressive)) != null && !usedKmers.contains(pk)) {
            stretchBuilder.insert(0, String.valueOf(pk.charAt(0)));

            tk = pk;
            usedKmers.add(pk);
        }

        return stretchBuilder.toString();
    }

    public static String getSeededStretchRight(CortexGraph cg, String kmer, int color, boolean beAggressive) {
        String tk = kmer;
        StringBuilder stretchBuilder = new StringBuilder();

        Set<String> usedKmers = new HashSet<>();

        String nk;
        while ((nk = CortexUtils.getNextKmer(cg, tk, color, beAggressive)) != null && !usedKmers.contains(nk)) {
            stretchBuilder.append(String.valueOf(nk.charAt(nk.length() - 1)));

            tk = nk;
            usedKmers.add(nk);
        }

        return stretchBuilder.toString();
    }

    public static String getSeededStretchLeft(CortexGraph clean, CortexGraph dirty, String kmer, int color, boolean beAggressive) {
        String tk = kmer;
        StringBuilder stretchBuilder = new StringBuilder();

        Set<String> usedKmers = new HashSet<>();

        String pk;
        while ((pk = CortexUtils.getPrevKmer(clean, dirty, tk, color, beAggressive)) != null && !usedKmers.contains(pk)) {
            stretchBuilder.insert(0, String.valueOf(pk.charAt(0)));

            tk = pk;
            usedKmers.add(pk);
        }

        return stretchBuilder.toString();
    }

    public static String getSeededStretchRight(CortexGraph clean, CortexGraph dirty, String kmer, int color, boolean beAggressive) {
        String tk = kmer;
        StringBuilder stretchBuilder = new StringBuilder();

        Set<String> usedKmers = new HashSet<>();

        String nk;
        while ((nk = CortexUtils.getNextKmer(clean, dirty, tk, color, beAggressive)) != null && !usedKmers.contains(nk)) {
            stretchBuilder.append(String.valueOf(nk.charAt(nk.length() - 1)));

            tk = nk;
            usedKmers.add(nk);
        }

        return stretchBuilder.toString();
    }

    /**
     * Given a kmer, walk left and right until we get to junctions on either end.
     *
     * @param cg    the multi-color graph
     * @param kmer  the start kmer
     * @param color the color to traverse
     * @return      a string containing the results of the walk
     */
    public static String getSeededStretch(CortexGraph cg, String kmer, int color) {
        return getSeededStretchLeft(cg, kmer, color, false) + kmer + getSeededStretchRight(cg, kmer, color, false);
    }

    public static String getSeededStretch(CortexGraph cg, String kmer, int color, boolean beAggressive) {
        return getSeededStretchLeft(cg, kmer, color, beAggressive) + kmer + getSeededStretchRight(cg, kmer, color, beAggressive);
    }

    public static String getSeededStretch(CortexGraph clean, CortexGraph dirty, String kmer, int color, boolean beAggressive) {
        return getSeededStretchLeft(clean, dirty, kmer, color, beAggressive) + kmer + getSeededStretchRight(clean, dirty, kmer, color, beAggressive);
    }

    public static boolean hasKmer(CortexGraph cg, String kmer, int color) {
        CortexRecord cr = cg.findRecord(new CortexKmer(kmer));

        return (cr != null && cr.getCoverage(color) > 0);
    }

    public static boolean isNovelKmer(CortexRecord cr, int color) {
        if (cr == null) {
            throw new IndianaException("Unable to test for novelty on a null Cortex record");
        }

        if (cr.getCoverage(color) == 0) {
            return false;
        }

        for (int c = 0; c < cr.getNumColors(); c++) {
            if (c != color && cr.getCoverage(c) != 0) {
                return false;
            }
        }

        return true;
    }

    public static boolean isNovelKmer(CortexGraph cg, String kmer, int color) {
        CortexRecord cr = cg.findRecord(new CortexKmer(kmer));

        if (cr == null) {
            throw new IndianaException("Unable to test for novelty on a kmer that's not in the graph: '" + kmer + "'");
        }

        if (cr.getCoverage(color) == 0) {
            return false;
        }

        for (int c = 0; c < cr.getNumColors(); c++) {
            if (c != color && cr.getCoverage(c) != 0) {
                return false;
            }
        }

        return true;
    }

    public static String getNovelStretch(CortexGraph cg, String kmer, int color, boolean beAggressive) {
        String tk = kmer;
        StringBuilder novelStretchBuilder = new StringBuilder(tk);

        Set<String> usedKmers = new HashSet<>();

        String pk;
        while ((pk = CortexUtils.getPrevKmer(cg, tk, color, beAggressive)) != null && isNovelKmer(cg, pk, color) && !usedKmers.contains(pk)) {
            novelStretchBuilder.insert(0, String.valueOf(pk.charAt(0)));

            tk = pk;
            usedKmers.add(pk);
        }

        tk = kmer;
        usedKmers.clear();

        String nk;
        while ((nk = CortexUtils.getNextKmer(cg, tk, color, beAggressive)) != null && isNovelKmer(cg, nk, color) && !usedKmers.contains(nk)) {
            novelStretchBuilder.append(String.valueOf(nk.charAt(nk.length() - 1)));

            tk = nk;
            usedKmers.add(nk);
        }

        return novelStretchBuilder.toString();
    }

    private static boolean isNovelKmer(CortexRecord cr) {
        int cov_c0 = cr.getCoverage(0);
        int cov_sum = 0;

        for (int c = 0; c < cr.getNumColors(); c++) {
            cov_sum += cr.getCoverage(c);
        }

        return cov_c0 > 0 && cov_c0 == cov_sum;
    }

    private static String formatAttributes(Map<String, Object> attrs) {
        List<String> attrArray = new ArrayList<>();

        for (String attr : attrs.keySet()) {
            String value = attrs.get(attr).toString();

            attrArray.add(attr + "=\"" + value + "\"");
        }

        return Joiner.on(" ").join(attrArray);
    }

    private static void printGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> g) {
        try {
            File f = new File("testgraph.dot");
            //File p = new File("testgraph.png");
            File p = new File("testgraph.pdf");

            PrintStream o = new PrintStream(f);

            String indent = "  ";

            o.println("digraph G {");

            for (AnnotatedVertex v : g.vertexSet()) {
                Map<String, Object> attrs = new TreeMap<>();
                attrs.put("label", "");
                attrs.put("fillcolor", v.isNovel() ? "red" : "white");
                attrs.put("style", "filled");

                o.println(indent + "\"" + v.getKmer() + "\" [ " + formatAttributes(attrs) + " ];");
            }

            String[] colors = new String[] { "red", "blue", "green" };

            for (AnnotatedEdge e : g.edgeSet()) {
                String s = g.getEdgeSource(e).getKmer();
                String t = g.getEdgeTarget(e).getKmer();

                for (int c = 0; c < 3; c++) {
                    if (e.isPresent(c)) {
                        Map<String, Object> attrs = new TreeMap<>();
                        attrs.put("color", colors[c]);

                        o.println(indent + "\"" + s + "\" -> \"" + t + "\" [ " + formatAttributes(attrs) + " ];");
                    }
                }
            }

            o.println("}");

            o.close();

            //System.out.println("dot -Tpng -o" + p.getAbsolutePath() + " " + f.getAbsolutePath());
            //Runtime.getRuntime().exec("dot -Tpng -o" + p.getAbsolutePath() + " " + f.getAbsolutePath());
            Runtime.getRuntime().exec("dot -Tpdf -o" + p.getAbsolutePath() + " " + f.getAbsolutePath());
        } catch (FileNotFoundException e) {
            throw new IndianaException("File not found", e);
        } catch (IOException e) {
            throw new IndianaException("IO exception", e);
        }
    }

    public static boolean hasExpectedOverlap(AnnotatedVertex av0, AnnotatedVertex av1) {
        String sk0 = av0.getKmer();
        String sk1 = av1.getKmer();

        return (sk0.substring(1, sk0.length()).equals(sk1.substring(0, sk0.length() - 1)));
    }

    private static TraversalStopper<AnnotatedVertex, AnnotatedEdge> instantiateStopper(Class<? extends TraversalStopper<AnnotatedVertex, AnnotatedEdge>> stopperClass) {
        try {
            return stopperClass.newInstance();
        } catch (InstantiationException e) {
            throw new IndianaException("Could not instantiate stopper: ", e);
        } catch (IllegalAccessException e) {
            throw new IndianaException("Illegal access while trying to instantiate stopper: ", e);
        }
    }

    public static int addVertexAndConnect(DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs, AnnotatedVertex cv, Map<Integer, Set<String>> prevKmers, Map<Integer, Set<String>> nextKmers) {
        int numVerticesBefore = dfs.vertexSet().size();

        dfs.addVertex(cv);

        Map<Integer, Set<String>> adjKmers;
        for (int c = 0; c < 3; c++) {
            for (int i = 0; i < 2; i++) {
                boolean edgesAreForward = i == 1;
                adjKmers = edgesAreForward ? nextKmers : prevKmers;

                for (String adjKmer : adjKmers.get(c)) {
                    AnnotatedVertex av = new AnnotatedVertex(adjKmer);

                    dfs.addVertex(av);

                    AnnotatedEdge aen;

                    if (edgesAreForward && dfs.containsEdge(cv, av)) {
                        aen = dfs.getEdge(cv, av);
                    } else if (!edgesAreForward && dfs.containsEdge(av, cv)) {
                        aen = dfs.getEdge(av, cv);
                    } else {
                        aen = new AnnotatedEdge();
                    }

                    aen.setPresence(c);

                    if (edgesAreForward) {
                        if (!hasExpectedOverlap(cv, av)) { throw new IndianaException("Kmers '" + cv.getKmer() + "' and '" + av.getKmer() + "' did not have expected overlap"); }

                        dfs.addEdge(cv, av, aen);
                    } else {
                        if (!hasExpectedOverlap(av, cv)) { throw new IndianaException("Kmers '" + cv.getKmer() + "' and '" + av.getKmer() + "' did not have expected overlap"); }

                        dfs.addEdge(av, cv, aen);
                    }
                }
            }
        }

        return dfs.vertexSet().size() - numVerticesBefore;
    }

    public static DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs(CortexGraph clean, CortexGraph dirty, String kmer, int color, DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, Class<? extends TraversalStopper<AnnotatedVertex, AnnotatedEdge>> stopperClass, int depth, boolean goForward, HashSet<String> history) {
        String firstKmer = new String(kmer);

        Map<Integer, Set<String>> sourceKmersAllColors = goForward ? CortexUtils.getPrevKmers(clean, dirty, kmer) : CortexUtils.getNextKmers(clean, dirty, kmer);
        Set<String> sourceKmers = sourceKmersAllColors.get(color);

        DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = new DefaultDirectedGraph<>(AnnotatedEdge.class);
        TraversalStopper<AnnotatedVertex, AnnotatedEdge> stopper = instantiateStopper(stopperClass);

        Map<Integer, Set<String>> adjKmers;

        do {
            AnnotatedVertex cv = new AnnotatedVertex(kmer);

            CortexRecord cr = clean.findRecord(new CortexKmer(cv.getKmer()));
            if (cr == null && dirty != null) {
                cr = dirty.findRecord(new CortexKmer(cv.getKmer()));
            }

            Map<Integer, Set<String>> prevKmers = CortexUtils.getPrevKmers(clean, dirty, cv.getKmer());
            Map<Integer, Set<String>> nextKmers = CortexUtils.getNextKmers(clean, dirty, cv.getKmer());
            adjKmers = goForward ? nextKmers : prevKmers;

            int numVerticesAdded = addVertexAndConnect(dfs, cv, prevKmers, nextKmers);

            if (stopper.keepGoing(cr, g, depth, dfs.vertexSet().size(), adjKmers.get(color).size()) && !sourceKmers.contains(kmer) && !history.contains(kmer)) {
                history.add(kmer);

                if (adjKmers.get(color).size() == 1) {
                    kmer = adjKmers.get(color).iterator().next();
                } else if (adjKmers.get(color).size() != 1) {
                    boolean childrenWereSuccessful = false;

                    for (String ak : adjKmers.get(color)) {
                        if (!ak.equals(firstKmer)) {
                            DirectedGraph<AnnotatedVertex, AnnotatedEdge> branch = dfs(clean, dirty, ak, color, g, stopperClass, depth + (CortexUtils.isNovelKmer(cr) ? 0 : 1), goForward, history);

                            if (branch != null) {
                                Graphs.addGraph(dfs, branch);
                                childrenWereSuccessful = true;
                            } else {
                                for (AnnotatedVertex av : dfs.vertexSet()) {
                                    if (av.getKmer().equals(ak)) {
                                        av.setFlag("branchRejected");
                                    }
                                }
                            }
                        }
                    }

                    if (childrenWereSuccessful || stopper.hasTraversalSucceeded(cr, g, depth, dfs.vertexSet().size(), 0)) {
                        return dfs;
                    } else {
                        for (AnnotatedVertex av : dfs.vertexSet()) {
                            if (av.getKmer().equals(kmer)) {
                                av.setFlag("branchRejected");
                            }
                        }
                    }
                }
            } else if (stopper.traversalSucceeded()) {
                return dfs;
            } else {
                return null;
            }
        } while (adjKmers.get(color).size() == 1);

        return null;
    }

    public static DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs(CortexGraph clean, CortexGraph dirty, String kmer, int color, DirectedGraph<AnnotatedVertex, AnnotatedEdge> sg0, Class<? extends TraversalStopper<AnnotatedVertex, AnnotatedEdge>> stopperClass) {
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = new DefaultDirectedGraph<>(AnnotatedEdge.class);

        DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfsf = dfs(clean, dirty, kmer, color, sg0, stopperClass, 0, true, new HashSet<>());
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfsb = dfs(clean, dirty, kmer, color, sg0, stopperClass, 0, false, new HashSet<>());

        if (dfsf != null) { Graphs.addGraph(dfs, dfsf); }
        if (dfsb != null) { Graphs.addGraph(dfs, dfsb); }

        return dfs;
    }

    /**
     * Given a kmer, extract the local subgraph by a depth-first search
     *
     * @param cg    the multi-color graph
     * @param kmer  the start kmer
     * @param color the color in which to begin the traversal
     * @return      a subgraph containing the local context
     */
    public static DirectedGraph<AnnotatedVertex, AnnotatedEdge> getSeededSubgraph(CortexGraph cg, String kmer, int color) {
        String stretch = CortexUtils.getSeededStretch(cg, kmer, color);

        DirectedGraph<AnnotatedVertex, AnnotatedEdge> sgc = new DefaultDirectedGraph<>(AnnotatedEdge.class);

        for (int i = 0; i <= stretch.length() - cg.getKmerSize() - 1; i++) {
            String sk0 = stretch.substring(i, i + cg.getKmerSize());
            String sk1 = stretch.substring(i + 1, i + 1 + cg.getKmerSize());

            AnnotatedVertex ak0 = new AnnotatedVertex(sk0);
            AnnotatedVertex ak1 = new AnnotatedVertex(sk1);

            CortexRecord cr0 = cg.findRecord(new CortexKmer(ak0.getKmer()));
            CortexRecord cr1 = cg.findRecord(new CortexKmer(ak1.getKmer()));

            if (isNovelKmer(cr0)) { ak0.setNovel(); }
            if (isNovelKmer(cr1)) { ak1.setNovel(); }

            sgc.addVertex(ak0);
            sgc.addVertex(ak1);
            sgc.addEdge(ak0, ak1, new AnnotatedEdge(true, false, false));
        }

        printGraph(sgc);

        Set<AnnotatedVertex> childVertices = new HashSet<>(sgc.vertexSet());

        for (int c = 1; c <= 2; c++) {
            for (AnnotatedVertex ak : childVertices) {
                CortexKmer ck = new CortexKmer(ak.getKmer());
                CortexRecord cr = cg.findRecord(ck);

                if (cr != null && cr.getCoverage(c) > 0) {
                    Set<String> nextKmers = CortexUtils.getNextKmers(cg, ak.getKmer(), c);

                    Stack<Pair<String, String>> kmerStack = new Stack<>();
                    for (String nextKmer : nextKmers) {
                        kmerStack.push(new Pair<>(ak.getKmer(), nextKmer));
                    }

                    while (!kmerStack.isEmpty()) {
                        Pair<String, String> p = kmerStack.pop();

                        AnnotatedVertex ak0 = new AnnotatedVertex(p.getFirst());
                        AnnotatedVertex ak1 = new AnnotatedVertex(p.getSecond());

                        if (!sgc.containsVertex(ak0)) {
                            sgc.addVertex(ak0);
                        }
                        if (!sgc.containsVertex(ak1)) {
                            sgc.addVertex(ak1);
                        }

                        if (!sgc.containsEdge(ak0, ak1)) {
                            sgc.addEdge(ak0, ak1, new AnnotatedEdge());
                        }

                        sgc.getEdge(ak0, ak1).set(c, true);

                        CortexRecord crNext = cg.findRecord(new CortexKmer(ak1.getKmer()));
                        nextKmers = CortexUtils.getNextKmers(cg, crNext.getKmerAsString(), c);

                        if (crNext.getCoverage(0) != 0) { // we've rejoined the child's graph!
                            for (String nextKmer : nextKmers) {
                                AnnotatedVertex akn = new AnnotatedVertex(nextKmer);

                                if (!sgc.containsVertex(akn)) { sgc.addVertex(akn); }

                                if (!sgc.containsEdge(ak1, akn)) {
                                    sgc.addEdge(ak1, akn, new AnnotatedEdge());
                                }

                                sgc.getEdge(ak1, akn).set(c, true);
                            }
                        } else {
                            for (String nextKmer : nextKmers) {
                                AnnotatedVertex akn = new AnnotatedVertex(nextKmer);

                                if (!sgc.containsVertex(akn)) {
                                    kmerStack.push(new Pair<>(p.getSecond(), nextKmer));
                                }
                            }
                        }

                        printGraph(sgc);
                    }
                }

                printGraph(sgc);

                if (cr != null && cr.getCoverage(c) > 0) {
                    Set<String> prevKmers = CortexUtils.getPrevKmers(cg, ak.getKmer(), c);

                    Stack<Pair<String, String>> kmerStack = new Stack<>();
                    for (String prevKmer : prevKmers) {
                        kmerStack.push(new Pair<>(prevKmer, ak.getKmer()));
                    }

                    while (!kmerStack.isEmpty()) {
                        Pair<String, String> p = kmerStack.pop();

                        AnnotatedVertex ak0 = new AnnotatedVertex(p.getFirst());
                        AnnotatedVertex ak1 = new AnnotatedVertex(p.getSecond());

                        if (!sgc.containsVertex(ak0)) { sgc.addVertex(ak0); }
                        if (!sgc.containsVertex(ak1)) { sgc.addVertex(ak1); }

                        if (!sgc.containsEdge(ak0, ak1)) {
                            sgc.addEdge(ak0, ak1, new AnnotatedEdge());
                        }

                        sgc.getEdge(ak0, ak1).set(c, true);

                        CortexRecord crPrev = cg.findRecord(new CortexKmer(ak0.getKmer()));
                        prevKmers = CortexUtils.getPrevKmers(cg, crPrev.getKmerAsString(), c);

                        if (crPrev.getCoverage(0) != 0) { // we've rejoined the child's graph!
                            for (String prevKmer : prevKmers) {
                                AnnotatedVertex akp = new AnnotatedVertex(prevKmer);

                                if (!sgc.containsVertex(akp)) { sgc.addVertex(akp); }

                                if (!sgc.containsEdge(akp, ak0)) {
                                    sgc.addEdge(akp, ak0, new AnnotatedEdge());
                                }

                                sgc.getEdge(akp, ak0).set(c, true);
                            }
                        } else {
                            for (String prevKmer : prevKmers) {
                                AnnotatedVertex akp = new AnnotatedVertex(prevKmer);

                                if (!sgc.containsVertex(akp)) {
                                    kmerStack.push(new Pair<>(prevKmer, p.getFirst()));
                                }
                            }
                        }

                        printGraph(sgc);
                    }
                }

                printGraph(sgc);
            }
        }

        return sgc;
    }

    /**
     * Flip the information in a links record to the opposite orientation
     *
     * @param clr  the Cortex links record
     * @return     the flipped Cortex links record
     */
    public static CortexLinksRecord flipLinksRecord(CortexLinksRecord clr) {
        String rsk = SequenceUtils.reverseComplement(clr.getKmerAsString());
        List<CortexJunctionsRecord> cjrs = new ArrayList<>();
        for (CortexJunctionsRecord cjr : clr.getJunctions()) {
            cjrs.add(flipJunctionsRecord(cjr));
        }

        return new CortexLinksRecord(rsk, cjrs);
    }

    /**
     * Flip the information in a junctions record to the opposite orientation
     *
     * @param cjr  the Cortex link junctions record
     * @return     the flipped Cortex link junctions record
     */
    public static CortexJunctionsRecord flipJunctionsRecord(CortexJunctionsRecord cjr) {
        return new CortexJunctionsRecord(!cjr.isForward(),
                                         cjr.getNumKmers(),
                                         cjr.getNumJunctions(),
                                         cjr.getCoverages(),
                                         SequenceUtils.complement(cjr.getJunctions()),
                                         SequenceUtils.reverseComplement(cjr.getSeq()),
                                         cjr.getJunctionPositions()
                   );
    }

    /**
     * Adjust the information in a links record to the opposite orientation
     *
     * @param sk   the kmer in desired orientation
     * @param clr  the Cortex link junctions record
     * @return     the reoriented Cortex links record
     */
    public static CortexLinksRecord orientLinksRecord(String sk, CortexLinksRecord clr) {
        List<CortexJunctionsRecord> cjrs = new ArrayList<>();
        for (CortexJunctionsRecord cjr : clr.getJunctions()) {
            cjrs.add(orientJunctionsRecord(sk, cjr));
        }

        return new CortexLinksRecord(sk, cjrs);
    }

    /**
     * Adjust the information in a junctions to facilitate navigation in contig orientation
     *
     * @param sk   the kmer in desired orientation
     * @param cjr  the Cortex link junctions record
     * @return     the reoriented Cortex link junctions record
     */
    public static CortexJunctionsRecord orientJunctionsRecord(String sk, CortexJunctionsRecord cjr) {
        CortexKmer ck = new CortexKmer(sk);
        boolean isForward = cjr.isForward();
        String junctions = cjr.getJunctions();

        if (!ck.isFlipped() && isForward) { // --> F
            // do nothing
        } else if (!ck.isFlipped() && !isForward) { // --> R
            junctions = SequenceUtils.complement(junctions);
        } else if (ck.isFlipped() && isForward) { // <-- F
            isForward = !isForward;
            junctions = SequenceUtils.complement(junctions);
        } else if (ck.isFlipped() && !isForward) { // <-- R
            isForward = !isForward;
        }

        return new CortexJunctionsRecord(isForward, cjr.getNumKmers(), cjr.getNumJunctions(), cjr.getCoverages(), junctions, cjr.getSeq(), cjr.getJunctionPositions());
    }

    /**
     * Return a list of all of the kmers in the graph between a starting kmer and the end of the junctions record
     * @param cg   the Cortex graph
     * @param sk   the kmer in desired orientation
     * @param cjr  the Cortex link junctions record
     * @return     the list of kmers
     */
    public static List<String> getKmersInLinkByNavigation(CortexGraph cg, String sk, CortexJunctionsRecord cjr) {
        cjr = orientJunctionsRecord(sk, cjr);

        int junctionsUsed = 0;
        String curKmer = sk;
        boolean isForward = cjr.isForward();
        String junctions = cjr.getJunctions();

        List<String> kmersInLink = new ArrayList<>();
        kmersInLink.add(sk);

        if (isForward) {
            Set<String> nextKmers = CortexUtils.getNextKmers(cg, sk, 0);

            while (nextKmers.size() > 0 && junctionsUsed < junctions.length()) {
                if (nextKmers.size() == 1) {
                    curKmer = nextKmers.iterator().next();
                    nextKmers = CortexUtils.getNextKmers(cg, curKmer, 0);

                    kmersInLink.add(curKmer);
                } else {
                    char nbase = junctions.charAt(junctionsUsed);

                    String expectedNextKmer = curKmer.substring(1, curKmer.length()) + nbase;

                    if (nextKmers.contains(expectedNextKmer)) {
                        curKmer = expectedNextKmer;
                        nextKmers = CortexUtils.getNextKmers(cg, curKmer, 0);

                        kmersInLink.add(expectedNextKmer);
                    } else {
                        throw new IndianaException("Junction record specified a navigation that conflicted with the graph: " + junctions + " " + nbase + " " + expectedNextKmer + " " + nextKmers);
                    }

                    junctionsUsed++;
                }
            }
        } else {
            Set<String> prevKmers = CortexUtils.getPrevKmers(cg, sk, 0);

            while (prevKmers.size() > 0 && junctionsUsed < junctions.length()) {
                if (prevKmers.size() == 1) {
                    curKmer = prevKmers.iterator().next();
                    prevKmers = CortexUtils.getPrevKmers(cg, curKmer, 0);

                    kmersInLink.add(0, curKmer);
                } else {
                    char pbase = junctions.charAt(junctionsUsed);

                    String expectedPrevKmer = pbase + curKmer.substring(0, curKmer.length() - 1);

                    if (prevKmers.contains(expectedPrevKmer)) {
                        curKmer = expectedPrevKmer;
                        prevKmers = CortexUtils.getPrevKmers(cg, curKmer, 0);

                        kmersInLink.add(0, expectedPrevKmer);
                    } else {
                        throw new IndianaException("Junction record specified a navigation that conflicted with the graph: " + junctions + " " + pbase + " " + expectedPrevKmer + " " + prevKmers);
                    }

                    junctionsUsed++;
                }
            }
        }

        return kmersInLink;
    }

    public static List<String> getKmersInLinkFromSeq(CortexGraph cg, String sk, CortexJunctionsRecord cjr) {
        List<String> expectedKmers = new ArrayList<>();
        String seq = cjr.getSeq();
        CortexKmer ck = new CortexKmer(sk);
        if ((!ck.isFlipped() && !cjr.isForward()) || (ck.isFlipped() && cjr.isForward())) {
            seq = SequenceUtils.reverseComplement(seq);
        }
        for (int i = 0; i <= seq.length() - cg.getKmerSize(); i++) {
            String kmer = seq.substring(i, i + cg.getKmerSize());
            expectedKmers.add(kmer);
        }

        return expectedKmers;
    }

    public static List<String> getKmersInLink(CortexGraph cg, String sk, CortexJunctionsRecord cjr) {
        return (cjr.getSeq() != null) ? getKmersInLinkFromSeq(cg, sk, cjr) : getKmersInLinkByNavigation(cg, sk, cjr);
    }

    public static Map<String, Integer> getKmersAndCoverageInLink(CortexGraph cg, String sk, CortexLinksRecord clr) {
        Map<String, Integer> kmersAndCoverage = new LinkedHashMap<>();

        for (CortexJunctionsRecord cjr : clr.getJunctions()) {
            List<String> kil = getKmersInLink(cg, sk, cjr);

            for (int i = 0; i < kil.size() - 1; i++) {
                String kila = kil.get(i);
                String kilb = kil.get(i+1);
                String id = "#" + kila + "_" + kilb;

                int cov = cjr.getCoverage(0);

                if (!kmersAndCoverage.containsKey(id)) {
                    kmersAndCoverage.put(id, cov);
                } else {
                    kmersAndCoverage.put(id, kmersAndCoverage.get(id) + cov);
                }
            }
        }

        return kmersAndCoverage;
    }

    public static DataFrame<String, String, Integer> getKmersAndCoverageHVStyle(CortexGraph cg, String sk, CortexLinksRecord clr) {
        DataFrame<String, String, Integer> hv = new DataFrame<>(0);

        for (CortexJunctionsRecord cjr : clr.getJunctions()) {
            List<String> kil = getKmersInLink(cg, sk, cjr);
            int cov = cjr.getCoverage(0);

            for (int i = 0; i < kil.size(); i++) {
                String kili = kil.get(i);

                for (int j = 0; j < kil.size(); j++) {
                    String kilj = kil.get(j);

                    if (i != j) {
                        hv.set(kili, kilj, hv.get(kili, kilj) + cov);
                    }
                }
            }
        }

        return hv;
    }

    public static byte binaryNucleotideToChar(long nucleotide) {
        switch ((int) nucleotide) {
            case 0: return 'A';
            case 1: return 'C';
            case 2: return 'G';
            case 3: return 'T';
            default:
                throw new RuntimeException("Nucleotide '" + nucleotide + "' is not a valid binary nucleotide");
        }
    }

    public static long charToBinaryNucleotide(byte b) {
        switch (b) {
            case 'A' : return 0;
            case 'C' : return 1;
            case 'G' : return 2;
            case 'T' : return 3;
            default:
                throw new RuntimeException("Nucleotide '" + b + "' is not a valid character nucleotide");
        }
    }

    public static byte[] decodeBinaryKmer(long[] kmer, int kmerSize, int kmerBits) {
        byte[] rawKmer = new byte[kmerSize];

        long[] binaryKmer = Arrays.copyOf(kmer, kmer.length);

        for (int i = 0; i < binaryKmer.length; i++) {
            binaryKmer[i] = reverse(binaryKmer[i]);
        }

        for (int i = kmerSize - 1; i >= 0; i--) {
            rawKmer[i] = binaryNucleotideToChar(binaryKmer[kmerBits - 1] & 0x3);

            shiftBinaryKmerByOneBase(binaryKmer, kmerBits);
        }

        return rawKmer;
    }

    private static int getKmerBits(byte[] kmer) {
        return (int) Math.ceil(((float) kmer.length)/32.0);
    }

    public static long[] encodeBinaryKmer(byte[] kmer, long[] oldbk) {
        int numBits = getKmerBits(kmer);
        long[] binaryKmer = new long[numBits];

        for (int b = 0; b < numBits; b++) {
            for (int i = kmer.length - 32*(b+1); i < kmer.length - 32*b; i++) {
                if (i >= 0) {
                    long nuc = charToBinaryNucleotide(kmer[i]);

                    binaryKmer[numBits - b - 1] |= nuc;
                }

                if (i < kmer.length - 32*b - 1) {
                    binaryKmer[numBits - b - 1] <<= 2;
                }
            }

            binaryKmer[numBits - b - 1] = reverse(binaryKmer[numBits - b - 1]);
        }

        return binaryKmer;
    }

    public static void shiftBinaryKmerByOneBase(long[] binaryKmer, int bitfields) {
        for(int i = bitfields - 1; i > 0; i--) {
            binaryKmer[i] >>>= 2;
            binaryKmer[i] |= (binaryKmer[i-1] << 62); // & 0x3
        }
        binaryKmer[0] >>>= 2;
    }

    public static long reverse(long x) {
        ByteBuffer bbuf = ByteBuffer.allocate(8);
        bbuf.order(ByteOrder.BIG_ENDIAN);
        bbuf.putLong(x);
        bbuf.order(ByteOrder.LITTLE_ENDIAN);

        return bbuf.getLong(0);
    }
}
