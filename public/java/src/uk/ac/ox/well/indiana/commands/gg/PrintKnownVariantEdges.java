package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import htsjdk.samtools.cram.structure.Container;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.assembly.cortex.CortexGraphWalker;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksMap;
import uk.ac.ox.well.indiana.utils.io.fasta.MultiIndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

public class PrintKnownVariantEdges extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference genome")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="bed", shortName="b", doc="Bed file")
    public File BED;

    @Argument(fullName="graph", shortName="g", doc="Graph (sorted)")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Graph links")
    public ArrayList<File> LINKS;

    @Argument(fullName="novelKmers", shortName="n", doc="Novel kmers")
    public FastaSequenceFile NOVEL_KMERS;

    @Output
    public File out;

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

    private void printGraph(DirectedGraph<String, Link> g, File f, String leftFlank, String rightFlank, Set<CortexKmer> novelKmers) {
        try {
            PrintStream o = new PrintStream(f);

            o.println("digraph G {");
            o.println("    graph [fontname = \"Courier\"];");
            o.println("    node [fontname = \"Courier\"];");
            o.println("    edge [fontname = \"Courier\"];");

            for (String v : g.vertexSet()) {
                CortexKmer ck = new CortexKmer(v);
                String fontColor = novelKmers.contains(ck) ? "red" : "black";

                if (v.equals(leftFlank) || v.equals(rightFlank)) {
                    o.println("    \"" + v + "\" [ style=\"filled\" fillcolor=\"red\" fontcolor=\"" + fontColor + "\" ];");
                } else {
                    o.println("    \"" + v + "\" [ fontcolor=\"" + fontColor + "\" ];");
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
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    private Set<String> callFromGraph(DirectedGraph<String, Link> g, String leftFlank, String rightFlank, Set<CortexKmer> novelKmers) {
        Set<CortexKmer> novelKmersInSubGraph = new HashSet<CortexKmer>();

        for (String kmer : g.vertexSet()) {
            CortexKmer ck = new CortexKmer(kmer);
            if (novelKmers.contains(ck)) {
                novelKmersInSubGraph.add(ck);
            }
        }

        for (CortexKmer ck : novelKmersInSubGraph) {
            String skFw = ck.getKmerAsString();
            String skRc = SequenceUtils.reverseComplement(skFw);

            String sk = null;
            if (g.containsVertex(skFw)) {
                sk = skFw;
            } else if (g.containsVertex(skRc)) {
                sk = skRc;
            }

            if (sk != null) {
                Set<String> alleles = new HashSet<String>();

                List<String> stretch = new ArrayList<String>();
                stretch.add(sk);

                String tk = sk;
                while (!tk.equals(leftFlank) && inDegreeOfColor(g, tk, 0) == 1) {
                    Link e = incomingEdgesOf(g, tk, 0).iterator().next();
                    tk = g.getEdgeSource(e);

                    stretch.add(0, tk);
                }

                tk = sk;
                while (!tk.equals(rightFlank) && outDegreeOfColor(g, tk, 0) == 1) {
                    Link e = outgoingEdgesOf(g, tk, 0).iterator().next();
                    tk = g.getEdgeTarget(e);

                    stretch.add(tk);
                }

                StringBuilder ab = new StringBuilder();
                for (String kmer : stretch) {
                    ab.append(kmer.charAt(kmer.length() - 1));
                }

                alleles.add(ab.toString());

                if (leftFlank.equals(stretch.get(0)) && rightFlank.equals(stretch.get(stretch.size() - 1))) {
                    List<String> kmerStretch = new ArrayList<String>();
                    kmerStretch.add(leftFlank);

                    String fk = leftFlank;
                    while (!fk.equals(rightFlank) && (outDegreeOfColor(g, fk, 1) == 1 || outDegreeOfColor(g, fk, 2) == 1)) {
                        Link e1 = outDegreeOfColor(g, fk, 1) == 1 ? outgoingEdgesOf(g, fk, 1).iterator().next() : null;
                        Link e2 = outDegreeOfColor(g, fk, 2) == 1 ? outgoingEdgesOf(g, fk, 2).iterator().next() : null;

                        Link e = (e1 != null) ? e1 : e2;

                        fk = g.getEdgeTarget(e);
                        kmerStretch.add(fk);
                    }

                    StringBuilder pb = new StringBuilder();
                    for (String kmer : kmerStretch) {
                        pb.append(kmer.charAt(kmer.length() - 1));
                    }

                    alleles.add(pb.toString());

                    return alleles;
                }
            }
        }

        return new HashSet<String>();

        /*
        for (CortexKmer ck : novelKmersInSubGraph) {
            if (!usedNovelKmers.contains(ck)) {
                String skFw = ck.getKmerAsString();
                String skRc = SequenceUtils.reverseComplement(skFw);

                String sk = null;
                if (g.containsVertex(skFw)) {
                    sk = skFw;
                } else if (g.containsVertex(skRc)) {
                    sk = skRc;
                }

                if (sk != null) {
                    List<String> novelKmerStretch = new ArrayList<String>();
                    novelKmerStretch.add(sk);

                    String pk = sk;
                    while (inDegreeOfColor(g, pk, 0) == 1 && (inDegreeOfColor(g, pk, 1) == 0 && inDegreeOfColor(g, pk, 2) == 0)) {
                        Link e = incomingEdgesOf(g, pk, 0).iterator().next();

                        pk = g.getEdgeSource(e);
                        novelKmerStretch.add(0, pk);

                        if (outDegreeOfColor(g, pk, 0) != 1) { break; }
                    }

                    String nk = sk;
                    while (outDegreeOfColor(g, nk, 0) == 1 && (outDegreeOfColor(g, nk, 1) == 0 && outDegreeOfColor(g, nk, 2) == 0)) {
                        Link e = outgoingEdgesOf(g, nk, 0).iterator().next();

                        nk = g.getEdgeTarget(e);
                        novelKmerStretch.add(nk);

                        if (inDegreeOfColor(g, nk, 0) != 1) { break; }
                    }

                    StringBuilder sb = new StringBuilder();
                    for (String kmer : novelKmerStretch) {
                        usedNovelKmers.add(new CortexKmer(kmer));

                        sb.append(kmer.charAt(kmer.length() - 1));
                    }

                    if (sb.length() > 0) {
                        alleles.add(sb.toString());
                    }

                    //String leftFlank = novelKmerStretch.get(0);
                    //String rightFlank = novelKmerStretch.get(novelKmerStretch.size() - 1);

                    for (int c = 1; c <= 2; c++) {
                        List<String> kmerStretch = new ArrayList<String>();
                        kmerStretch.add(leftFlank);

                        String tk = leftFlank;
                        while (!tk.equals(rightFlank) && outDegreeOfColor(g, tk, c) == 1) {
                            Link e = outgoingEdgesOf(g, tk, c).iterator().next();

                            tk = g.getEdgeTarget(e);
                            kmerStretch.add(tk);
                        }

                        StringBuilder pb = new StringBuilder();
                        for (String kmer : kmerStretch) {
                            pb.append(kmer.charAt(kmer.length() - 1));
                        }

                        alleles.add(pb.toString());
                    }
                }
            }
        }
        */
    }

    @Override
    public void execute() {
        Set<CortexKmer> novelKmers = loadNovelKmers();
        int kmerSize = novelKmers.iterator().next().length();

        Map<String, CortexLinksMap> linksMap = new HashMap<String, CortexLinksMap>();
        for (File ctp : LINKS) {
            CortexLinksMap clm = new CortexLinksMap(ctp);
            linksMap.put(clm.getCortexLinks().getColor(0).getSampleName(), clm);
        }

        TableReader tr = new TableReader(BED, "chr", "start", "stop", "info");

        Map<String, Integer> eventTypes = new HashMap<String, Integer>();
        int count = 0, full = 0;
        for (Map<String, String> te : tr) {
            String chr = te.get("chr");
            int start = Integer.valueOf(te.get("start")) + 1;
            int stop = Integer.valueOf(te.get("stop"));
            Map<String, String> infoMap = parseInfoField(te.get("info"));

            String eventType = infoMap.get("denovo");
            ContainerUtils.increment(eventTypes, eventType);

            if (!eventType.equals("unknown")) {
                String seq = new String(REF.getSubsequenceAt(chr, start, stop).getBases());
                String alt = infoMap.get("alt");

                // Display
                int window = 100;

                int contextStart = (start - window >= 1) ? start - window : 1;
                int contextEnd = stop + window < REF.getSequence(chr).length() ? stop + window : REF.getSequence(chr).length() - 1;
                String context = new String(REF.getSubsequenceAt(chr, contextStart, contextEnd).getBases());
                StringBuilder novelty = new StringBuilder();
                StringBuilder edgesOut = new StringBuilder(StringUtil.repeatCharNTimes(' ', contextEnd - contextStart + kmerSize));
                StringBuilder edgesIn = new StringBuilder(StringUtil.repeatCharNTimes(' ', contextEnd - contextStart + kmerSize));

                String leftFlank = new String(REF.getSubsequenceAt(chr, start - kmerSize, start - 1).getBases());
                String rightFlank = new String(REF.getSubsequenceAt(chr, stop + 1, stop + kmerSize).getBases());

                //String leftFlank = new String(REF.getSubsequenceAt(chr, start - 1, start - 1 + kmerSize - 1).getBases());
                //String rightFlank = new String(REF.getSubsequenceAt(chr, stop + 1 - kmerSize, stop).getBases());

                for (int i = 0; i <= context.length() - kmerSize; i++) {
                    CortexKmer kmer = new CortexKmer(context.substring(i, i + kmerSize));
                    if (novelKmers.contains(kmer)) {
                        novelty.append(".");
                    } else {
                        novelty.append(" ");
                    }

                    CortexRecord cr = GRAPH.findRecord(kmer);
                    if (!kmer.isFlipped()) {
                        if (cr.getOutEdgesAsBytes(0).size() > 1) { edgesOut.setCharAt(i + kmerSize, '/'); }
                        else { edgesOut.setCharAt(i + kmerSize, ' '); }

                        if (cr.getInEdgesAsBytes(0).size() > 1) { edgesIn.setCharAt(i, '/'); }
                        else { edgesIn.setCharAt(i, ' '); }
                    } else {
                        if (cr.getOutEdgesAsBytes(0).size() > 1) { edgesIn.setCharAt(i, '/'); }
                        else { edgesIn.setCharAt(i, ' '); }

                        if (cr.getInEdgesAsBytes(0).size() > 1) { edgesOut.setCharAt(i + kmerSize, '/'); }
                        else { edgesOut.setCharAt(i + kmerSize, ' '); }
                    }
                }

                DirectedGraph<String, Link> g = new DefaultDirectedGraph<String, Link>(Link.class);

                // Build child's sequence context
                StringBuilder childSequenceLeft = null;
                Set<String> thisKmer = new HashSet<String>(Arrays.asList(leftFlank));
                String lk = null;
                while (thisKmer.size() == 1) {
                    String tk = thisKmer.iterator().next();

                    if (childSequenceLeft == null) { childSequenceLeft = new StringBuilder(tk); }
                    else { childSequenceLeft.insert(0, tk.charAt(0)); }

                    g.addVertex(tk);
                    if (lk != null) { g.addEdge(tk, lk, new Link()); }
                    lk = tk;

                    thisKmer = CortexUtils.getPrevKmers(GRAPH, tk, 0);
                }
                for (String tk1 : thisKmer) {
                    g.addVertex(tk1);
                    if (lk != null) { g.addEdge(tk1, lk, new Link()); }
                }
                //printGraph(g, out, leftFlank, rightFlank, novelKmers);

                thisKmer = new HashSet<String>(Arrays.asList(leftFlank));
                lk = null;
                while (thisKmer.size() == 1) {
                    String tk = thisKmer.iterator().next();

                    if (childSequenceLeft == null) { childSequenceLeft = new StringBuilder(tk); }
                    else { childSequenceLeft.append(tk.charAt(tk.length() - 1)); }

                    g.addVertex(tk);
                    if (lk != null) { g.addEdge(lk, tk, new Link()); }
                    lk = tk;

                    thisKmer = CortexUtils.getNextKmers(GRAPH, tk, 0);
                }
                for (String tk1 : thisKmer) {
                    g.addVertex(tk1);
                    if (lk != null) { g.addEdge(lk, tk1, new Link()); }
                }
                //printGraph(g, out, leftFlank, rightFlank, novelKmers);

                StringBuilder childSequenceRight = null;
                thisKmer = new HashSet<String>(Arrays.asList(rightFlank));
                lk = null;
                while (thisKmer.size() == 1) {
                    String tk = thisKmer.iterator().next();

                    if (childSequenceRight == null) { childSequenceRight = new StringBuilder(tk); }
                    else { childSequenceRight.insert(0, tk.charAt(0)); }

                    g.addVertex(tk);
                    if (lk != null) { g.addEdge(tk, lk, new Link()); }
                    lk = tk;

                    thisKmer = CortexUtils.getPrevKmers(GRAPH, tk, 0);
                }
                for (String tk1 : thisKmer) {
                    g.addVertex(tk1);
                    if (lk != null) { g.addEdge(tk1, lk, new Link()); }
                }
                //printGraph(g, out, leftFlank, rightFlank, novelKmers);

                thisKmer = new HashSet<String>(Arrays.asList(rightFlank));
                lk = null;
                while (thisKmer.size() == 1) {
                    String tk = thisKmer.iterator().next();

                    if (childSequenceRight == null) { childSequenceRight = new StringBuilder(tk); }
                    else { childSequenceRight.append(tk.charAt(tk.length() - 1)); }

                    g.addVertex(tk);
                    if (lk != null) { g.addEdge(lk, tk, new Link()); }
                    lk = tk;

                    thisKmer = CortexUtils.getNextKmers(GRAPH, tk, 0);
                }
                for (String tk1 : thisKmer) {
                    g.addVertex(tk1);
                    if (lk != null) { g.addEdge(lk, tk1, new Link()); }
                }

                printGraph(g, out, leftFlank, rightFlank, novelKmers);

                for (String gk : g.vertexSet()) {
                    if (gk.length() != kmerSize) {
                        log.info("1 {}: {}", gk, gk.length());
                    }
                }

                for (int c = 0; c <= 2; c++) {
                    log.debug("event: {}",   eventType);
                    log.debug("  seq: {}{}", StringUtil.repeatCharNTimes(' ', start - contextStart), seq);
                    log.debug("  alt: {}{}", StringUtil.repeatCharNTimes(' ', start - contextStart), alt);
                    log.debug("  out: {}{}", StringUtil.repeatCharNTimes(' ', kmerSize), edgesOut);
                    log.debug("  con: {}",   context);
                    log.debug("   in: {}",   edgesIn);
                    log.debug("  nov: {}",   novelty);
                    log.debug("   lf: {}{}", StringUtil.repeatCharNTimes(' ', start - kmerSize - contextStart), leftFlank);
                    log.debug("   rf: {}{}", StringUtil.repeatCharNTimes(' ', stop + 1 - contextStart), rightFlank);
                    //log.debug("   lf: {}{}", StringUtil.repeatCharNTimes(' ', start - 1 - contextStart), leftFlank);
                    //log.debug("   rf: {}{}", StringUtil.repeatCharNTimes(' ', stop + 1 - kmerSize - contextStart), rightFlank);
                    log.debug("     : {}",   childSequenceLeft);
                    log.debug("     : {}",   childSequenceRight);

                    Set<String> vertices = new HashSet<String>(g.vertexSet());

                    for (String sk : vertices) {
                        CortexKmer ck = new CortexKmer(sk);
                        String sampleName = GRAPH.getColor(c).getSampleName();

                        if (linksMap.get(sampleName).containsKey(ck)) {
                            for (CortexJunctionsRecord cjr : linksMap.get(sampleName).get(ck).getJunctions()) {
                                String linkSeq = cjr.getSeq();

                                if (ck.isFlipped() != cjr.isFlipped()) {
                                    linkSeq = SequenceUtils.reverseComplement(linkSeq);
                                }

                                //if (linkSeq.contains(leftFlank) || linkSeq.contains(rightFlank)) {
                                if (linkSeq.contains(alt)) {
                                    Map<String, Integer> kmerTotalCount = new HashMap<String, Integer>();

                                    int[] indices = new int[linkSeq.length() - kmerSize + 1];
                                    for (int q = 0; q <= linkSeq.length() - kmerSize; q++) {
                                        String tk = linkSeq.substring(q, q + kmerSize);

                                        ContainerUtils.increment(kmerTotalCount, tk);
                                        indices[q] = kmerTotalCount.get(tk);
                                    }

                                    for (int q = 0; q < linkSeq.length() - kmerSize; q++) {
                                        String tk = linkSeq.substring(q, q + kmerSize);
                                        String nk = linkSeq.substring(q + 1, q + 1 + kmerSize);

                                        int tkIndex = indices[q];
                                        int nkIndex = indices[q + 1];

                                        String tki = tk + (tkIndex == 1 ? "" : "." + tkIndex);
                                        String nki = nk + (nkIndex == 1 ? "" : "." + nkIndex);

                                        g.addVertex(tki);
                                        g.addVertex(nki);

                                        if (!g.containsEdge(tki, nki)) {
                                            g.addEdge(tki, nki, new Link(c, cjr.getCoverage(0)));
                                        } else {
                                            g.getEdge(tki, nki).addCoverage(c, cjr.getCoverage(0));
                                        }

                                        //printGraph(g, out, leftFlank, rightFlank, novelKmers);
                                        //log.info("");
                                    }
                                }
                            }
                        }
                    }
                }

                for (String gk : g.vertexSet()) {
                    if (gk.length() != kmerSize) {
                        log.info("2 {}: {}", gk, gk.length());
                    }
                }

                String foundLeftFlank = context.substring(100, 100 + kmerSize);
                while (inDegreeOfColor(g, foundLeftFlank, 0) == 1 && (outDegreeOfColor(g, foundLeftFlank, 1) == 0 || outDegreeOfColor(g, foundLeftFlank, 2) == 0)) {
                    Link e = incomingEdgesOf(g, foundLeftFlank, 0).iterator().next();
                    foundLeftFlank = g.getEdgeSource(e);
                }

                String foundRightFlank = context.substring(context.length() - 101, context.length() - 101 + kmerSize);
                while (outDegreeOfColor(g, foundRightFlank, 0) == 1 && (inDegreeOfColor(g, foundRightFlank, 1) == 0 || inDegreeOfColor(g, foundRightFlank, 2) == 0)) {
                    Link e = outgoingEdgesOf(g, foundRightFlank, 0).iterator().next();
                    foundRightFlank = g.getEdgeTarget(e);
                }

                printGraph(g, out, foundLeftFlank, foundRightFlank, novelKmers);
                Set<String> alleles = callFromGraph(g, foundLeftFlank, foundRightFlank, novelKmers);

                log.info("Allele: {}", Joiner.on(", ").join(alleles));
                log.info("{} {} {}:{}-{} {}", count, full, chr, start, stop, eventType);

                //break;
            }
        }

        log.info("  {}", Joiner.on("; ").withKeyValueSeparator("=").join(eventTypes));
    }
}
