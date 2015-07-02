package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.EdgeReversedGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public class PrintNovelSubgraph extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="graphRaw", shortName="gr", doc="Graph (raw)", required=false)
    public CortexGraph GRAPH_RAW;

    @Argument(fullName="novelGraph", shortName="n", doc="Graph of novel kmers")
    public CortexGraph NOVEL;

    @Output
    public File out;

    private String formatAttributes(Map<String, Object> attrs) {
        List<String> attrArray = new ArrayList<String>();

        for (String attr : attrs.keySet()) {
            String value = attrs.get(attr).toString();

            attrArray.add(attr + "=\"" + value + "\"");
        }

        return Joiner.on(" ").join(attrArray);
    }

    public void printGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> g, int stretchNum) {
        try {
            File f = new File(out.getAbsolutePath() + ".stretch" + stretchNum + ".dot");
            File p = new File(out.getAbsolutePath() + ".stretch" + stretchNum + ".png");

            PrintStream o = new PrintStream(f);

            String indent = "  ";

            o.println("digraph G {");

            for (AnnotatedVertex v : g.vertexSet()) {
                Map<String, Object> attrs = new TreeMap<String, Object>();
                //attrs.put("label", "");
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
                        Map<String, Object> attrs = new TreeMap<String, Object>();
                        attrs.put("color", colors[c]);

                        o.println(indent + "\"" + s + "\" -> \"" + t + "\" [ " + formatAttributes(attrs) + " ];");
                    }
                }
            }

            o.println("}");

            o.close();

            log.info("dot -Tpng -o" + p.getAbsolutePath() + " " + f.getAbsolutePath());
            Runtime.getRuntime().exec("dot -Tpng -o" + p.getAbsolutePath() + " " + f.getAbsolutePath());
        } catch (FileNotFoundException e) {
            throw new IndianaException("File not found", e);
        } catch (IOException e) {
            throw new IndianaException("IO exception", e);
        }
    }

    private void addGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, DirectedGraph<String, DefaultEdge> g, int color, Map<CortexKmer, Boolean> novelKmers) {
        for (String v : g.vertexSet()) {
            AnnotatedVertex av = new AnnotatedVertex(v);
            if (novelKmers.containsKey(new CortexKmer(v))) {
                av.setNovel();
            }

            a.addVertex(av);
        }

        for (DefaultEdge e : g.edgeSet()) {
            String s0 = g.getEdgeSource(e);
            String s1 = g.getEdgeTarget(e);

            AnnotatedVertex a0 = new AnnotatedVertex(s0);
            AnnotatedVertex a1 = new AnnotatedVertex(s1);

            if (!a.containsEdge(a0, a1)) {
                a.addEdge(a0, a1, new AnnotatedEdge());
            }

            a.getEdge(a0, a1).set(color, true);
        }
    }

    private Set<AnnotatedVertex> getListOfVerticesToRemove(DirectedGraph<AnnotatedVertex, AnnotatedEdge> a, DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> dfs) {
        Set<AnnotatedVertex> toRemove = new HashSet<AnnotatedVertex>();

        while (dfs.hasNext()) {
            AnnotatedVertex v1 = dfs.next();

            if (!v1.isNovel() && !toRemove.contains(v1)) {
                Set<AnnotatedEdge> oes = a.outgoingEdgesOf(v1);

                for (AnnotatedEdge oe : oes) {
                    if (oe.isPresent(0) && (oe.isPresent(1) || oe.isPresent(2))) {
                        AnnotatedVertex vt = a.getEdgeTarget(oe);

                        DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> subdfs = new DepthFirstIterator<AnnotatedVertex, AnnotatedEdge>(a, vt);

                        while (subdfs.hasNext()) {
                            AnnotatedVertex vs = subdfs.next();

                            toRemove.add(vs);
                        }
                    }
                }
            }
        }

        return toRemove;
    }

    private void trimGraph(DirectedGraph<AnnotatedVertex, AnnotatedEdge> af) {
        DirectedGraph<AnnotatedVertex, AnnotatedEdge> ar = new EdgeReversedGraph<AnnotatedVertex, AnnotatedEdge>(af);

        Set<AnnotatedVertex> toRemove = new HashSet<AnnotatedVertex>();

        for (AnnotatedVertex av : af.vertexSet()) {
            if (av.isNovel()) {
                AnnotatedVertex v0 = av;

                DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> dfsf = new DepthFirstIterator<AnnotatedVertex, AnnotatedEdge>(af, v0);
                DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> dfsr = new DepthFirstIterator<AnnotatedVertex, AnnotatedEdge>(ar, v0);

                Set<AnnotatedVertex> toRemoveFwd = getListOfVerticesToRemove(af, dfsf);
                Set<AnnotatedVertex> toRemoveRev = getListOfVerticesToRemove(ar, dfsr);

                log.info("  fwd: {}/{}", toRemoveFwd.size(), af.vertexSet().size());
                log.info("  rev: {}/{}", toRemoveRev.size(), af.vertexSet().size());

                toRemove.addAll(toRemoveFwd);
                toRemove.addAll(toRemoveRev);
            }
        }

        af.removeAllVertices(toRemove);
    }

    @Override
    public void execute() {
        Map<CortexKmer, Boolean> novelKmers = new HashMap<CortexKmer, Boolean>();
        for (CortexRecord cr : NOVEL) {
            novelKmers.put(new CortexKmer(cr.getKmerAsString()), true);
        }

        int totalNovelKmersUsed = 0;
        int stretchNum = 0;
        for (CortexKmer novelKmer : novelKmers.keySet()) {
            if (novelKmers.get(novelKmer)) {
                int novelKmersUsed = 0;

                // do stuff
                DirectedGraph<String, DefaultEdge> sg0 = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
                DirectedGraph<String, DefaultEdge> sg1 = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
                DirectedGraph<String, DefaultEdge> sg2 = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);

                String stretch = CortexUtils.getSeededStretch(GRAPH, novelKmer.getKmerAsString(), 0);

                log.info("  stretch: {}", stretchNum);
                log.info("    processing child graph...");

                for (int i = 0; i <= stretch.length() - novelKmer.length(); i++) {
                    String kmer = stretch.substring(i, i + novelKmer.length());

                    Graphs.addGraph(sg0, CortexUtils.dfs(GRAPH, kmer, 0, null, new AbstractTraversalStopper() {
                        @Override
                        public boolean hasTraversalSucceeded(CortexRecord cr, DirectedGraph<String, DefaultEdge> g, int depth) {
                            //return cr.getCoverage(1) != 0 || cr.getCoverage(2) != 0;
                            return false;
                        }

                        @Override
                        public boolean hasTraversalFailed(CortexRecord cr, DirectedGraph<String, DefaultEdge> g, int junctions) {
                            return junctions >= maxJunctions;
                        }
                    }));
                }

                //
                for (int c = 1; c <= 2; c++) {
                    log.info("    processing parent {} graph... ({} kmers)", c, sg0.vertexSet().size());

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
                                    return junctions >= maxJunctions;
                                }
                            };

                            Graphs.addGraph(sg, CortexUtils.dfs(GRAPH, kmer, c, sg0, stopper));
                        }
                    }
                }
                //

                log.info("    combining graphs...");

                DirectedGraph<AnnotatedVertex, AnnotatedEdge> ag = new DefaultDirectedGraph<AnnotatedVertex, AnnotatedEdge>(AnnotatedEdge.class);

                addGraph(ag, sg0, 0, novelKmers);
                addGraph(ag, sg1, 1, novelKmers);
                addGraph(ag, sg2, 2, novelKmers);

                log.info("    trimming graph...");
                trimGraph(ag);

                printGraph(ag, 100);
                printGraph(ag, stretchNum);

                for (AnnotatedVertex ak : ag.vertexSet()) {
                    CortexKmer ck = new CortexKmer(ak.getKmer());

                    if (novelKmers.containsKey(ck) && novelKmers.get(ck)) {
                        totalNovelKmersUsed++;
                        novelKmersUsed++;
                        novelKmers.put(ck, false);
                    }
                }

                stretchNum++;
            }
        }

        log.info("Num stretches: {}", stretchNum);
    }
}
