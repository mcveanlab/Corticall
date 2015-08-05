package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.ProcessExecutor;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

public class GraphGenotyperOld extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph (sorted)")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Graph links")
    public ArrayList<File> LINKS;

    @Argument(fullName="novelKmers", shortName="n", doc="Novel kmers")
    public FastaSequenceFile NOVEL_KMERS;

    @Output
    public PrintStream out;

    private class Link extends DefaultEdge {
        public Map<Integer, Integer> covMap = new HashMap<Integer, Integer>();

        public Link() {
            covMap.put(-1, 1);
        }

        public Link(int color, int coverage) {
            addCoverage(color, coverage);
        }

        public void addCoverage(int color, int coverage) {
            if (!covMap.containsKey(color)) {
                covMap.put(color, 0);
            }

            covMap.put(color, covMap.get(color) + coverage);
        }

        public boolean hasEdgeInColor(int color) {
            return (covMap.containsKey(color) && covMap.get(color) > 0);
        }
    }

    private int inDegreeOfColor(DirectedGraph<String, Link> g, String kmer, int color) {
        int inDegree = 0;

        for (Link e : g.incomingEdgesOf(kmer)) {
            if (e.hasEdgeInColor(color)) { inDegree++; }
        }

        return inDegree;
    }

    private int outDegreeOfColor(DirectedGraph<String, Link> g, String kmer, int color) {
        int outDegree = 0;

        for (Link e : g.outgoingEdgesOf(kmer)) {
            if (e.hasEdgeInColor(color)) { outDegree++; }
        }

        return outDegree;
    }

    private Set<Link> incomingEdgesOf(DirectedGraph<String, Link> g, String kmer, int color) {
        Set<Link> edges = new HashSet<Link>();

        for (Link e : g.incomingEdgesOf(kmer)) {
            if (e.hasEdgeInColor(color)) {
                edges.add(e);
            }
        }

        return edges;
    }

    private Set<Link> outgoingEdgesOf(DirectedGraph<String, Link> g, String kmer, int color) {
        Set<Link> edges = new HashSet<Link>();

        for (Link e : g.outgoingEdgesOf(kmer)) {
            if (e.hasEdgeInColor(color)) {
                edges.add(e);
            }
        }

        return edges;
    }

    private Map<String, String> parseInfoField(String info) {
        Map<String, String> infoMap = new HashMap<String, String>();
        String[] keyValuePairs = info.split(";");

        for (String keyValuePair : keyValuePairs) {
            String[] kv = keyValuePair.split("=");
            infoMap.put(kv[0], kv[1]);
        }

        return infoMap;
    }

    private Set<CortexKmer> loadNovelKmers() {
        Set<CortexKmer> novelKmers = new HashSet<CortexKmer>();
        ReferenceSequence rseq;
        while ((rseq = NOVEL_KMERS.nextSequence()) != null) {
            CortexKmer ck = new CortexKmer(rseq.getBases());

            novelKmers.add(ck);
        }

        return novelKmers;
    }

    private void printGraph(DirectedGraph<String, Link> g, File f, Set<String> leftFlankKmers, Set<String> rightFlankKmers, Set<CortexKmer> variantKmers, Set<CortexKmer> novelKmers) {
        try {
            PrintStream o = new PrintStream(f);

            o.println("digraph G {");
            o.println("    graph [fontname = \"Courier\"];");
            o.println("    node [shape=point label=\"\" fontname = \"Courier\"];");
            o.println("    edge [fontname = \"Courier\"];");

            for (String v : g.vertexSet()) {
                CortexKmer ck = new CortexKmer(v);

                if (novelKmers.contains(ck)) {
                    o.println("    \"" + v + "\" [ label=\"" + v + "\" shape=rect style=\"filled\" fillcolor=\"red\" ];");
                } else if (variantKmers.contains(ck)) {
                    o.println("    \"" + v + "\" [ label=\"" + v + "\" shape=rect style=\"filled\" fillcolor=\"cyan\" ];");
                } else if (leftFlankKmers.contains(v) || rightFlankKmers.contains(v)) {
                    o.println("    \"" + v + "\" [ label=\"" + v + "\" shape=rect style=\"filled\" fillcolor=\"orange\" ];");
                } else {
                    o.println("    \"" + v + "\" [ label=\"" + v.charAt(v.length() - 1) + "\" shape=rect ];");
                }
            }

            Map<Integer, String> colorMap = new HashMap<Integer, String>();
            colorMap.put(-1, "black");
            colorMap.put(0,  "red");
            colorMap.put(1,  "green");
            colorMap.put(2,  "blue");

            for (Link e : g.edgeSet()) {
                String s = g.getEdgeSource(e);
                String t = g.getEdgeTarget(e);

                Link l = g.getEdge(s, t);

                for (Integer color : l.covMap.keySet()) {
                    int cov = l.covMap.get(color);
                    String col = colorMap.get(color);
                    //o.println("    \"" + s + "\" -> \"" + t + "\" [ penwidth=" + Math.log10(cov + 1) + " color=\"" + col + "\" label=\"" + cov + "\" ];");
                    o.println("    \"" + s + "\" -> \"" + t + "\" [ penwidth=" + Math.log10(cov + 1) + " color=\"" + col + "\" ];");
                }
            }

            o.println("}");

            o.close();

            ProcessExecutor.execute("dot -Tpdf -otestgraph.pdf testgraph.dot");
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void execute() {
        Set<CortexKmer> novelKmers = loadNovelKmers();

        // Start by putting novel kmers into the longest stretches we can make
        log.info("Finding linear stretches of novel sequence...");
        List<List<String>> stretches = getStretchesOfNovelKmers(novelKmers);

        // For each stretch, look for the parental branch points nearby
        for (int c = 1; c <= 2; c++) {
            for (List<String> stretch : stretches) {
                String prevKmer = stretch.get(0);

                do {
                    prevKmer = CortexUtils.getPrevKmer(GRAPH, prevKmer);
                } while (prevKmer != null && CortexUtils.getNextKmer(GRAPH, prevKmer, c, false) == null);

                String nextKmer = stretch.get(stretch.size() - 1);

                do {
                    nextKmer = CortexUtils.getNextKmer(GRAPH, nextKmer);
                } while (nextKmer != null && CortexUtils.getPrevKmer(GRAPH, nextKmer, c, false) == null);

                String prevKmerToNextChild  = CortexUtils.getNextKmer(GRAPH, prevKmer, 0, false);
                String prevKmerToNextParent = CortexUtils.getNextKmer(GRAPH, prevKmer, c, false);

                String nextKmerToPrevChild  = CortexUtils.getPrevKmer(GRAPH, nextKmer, 0, false);
                String nextKmerToPrevParent = CortexUtils.getPrevKmer(GRAPH, nextKmer, c, false);

                Map<Integer, String> alleles = new TreeMap<Integer, String>();

                if (prevKmerToNextChild != null && prevKmerToNextParent != null && !prevKmerToNextChild.equals(prevKmerToNextParent) &&
                    nextKmerToPrevChild != null && nextKmerToPrevParent != null && !nextKmerToPrevChild.equals(nextKmerToPrevParent)
                ) {
                    for (Integer ca : Arrays.asList(0, c)) {
                        String tk = prevKmer;
                        StringBuilder sf = new StringBuilder(tk);
                        boolean fComplete = false;

                        while ((tk = CortexUtils.getNextKmer(GRAPH, tk, ca, false)) != null) {
                            if (tk.equals(nextKmer)) {
                                fComplete = true;
                                break;
                            }

                            sf.append(tk.charAt(tk.length() - 1));
                        }

                        tk = nextKmer;
                        StringBuilder sr = new StringBuilder(tk);
                        boolean rComplete = false;

                        while ((tk = CortexUtils.getPrevKmer(GRAPH, tk, ca, false)) != null) {
                            if (tk.equals(prevKmer)) {
                                rComplete = true;
                                break;
                            }

                            sr.append(tk.charAt(0));
                        }

                        if      (fComplete) { alleles.put(ca, sf.toString()); }
                        else if (rComplete) { alleles.put(ca, sr.toString()); }
                    }
                }

                log.info("  stretch: {}, alleles: {}", stretch.size(), Joiner.on(", ").withKeyValueSeparator("=").join(alleles));
            }
        }

        log.info("Novel kmers: {}", novelKmers.size());
        log.info("  Stretches: {}", stretches.size());
        log.info("");
    }

    private List<List<String>> getStretchesOfNovelKmers(Set<CortexKmer> novelKmers) {
        List<List<String>> stretches = new ArrayList<List<String>>();
        Set<CortexKmer> usedNovelKmers = new HashSet<CortexKmer>();

        for (CortexKmer ck : novelKmers) {
            if (!usedNovelKmers.contains(ck)) {
                String sk = ck.getKmerAsString();

                List<String> stretch = new ArrayList<String>(Arrays.asList(sk));

                String tk = sk;
                while ((tk = CortexUtils.getPrevKmer(GRAPH, tk, 0, false)) != null && novelKmers.contains(new CortexKmer(tk))) {
                    stretch.add(0, tk);
                }

                tk = sk;
                while ((tk = CortexUtils.getNextKmer(GRAPH, tk, 0, false)) != null && novelKmers.contains(new CortexKmer(tk))) {
                    stretch.add(tk);
                }

                stretches.add(stretch);

                for (String kmer : stretch) {
                    CortexKmer ck0 = new CortexKmer(kmer);

                    if (novelKmers.contains(ck0)) {
                        usedNovelKmers.add(ck0);
                    }
                }

                if (stretches.size() % 1000 == 0) {
                    log.info("  stretches: {}, novel kmers: {}/{}", stretches.size(), usedNovelKmers.size(), novelKmers.size());
                }
            }
        }

        return stretches;
    }
}
