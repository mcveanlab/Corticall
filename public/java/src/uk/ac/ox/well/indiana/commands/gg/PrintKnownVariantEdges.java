package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
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

        public void addCoverage(int color, int coverage) {
            if (!covMap.containsKey(color)) {
                covMap.put(color, 0);
            }

            covMap.put(color, covMap.get(color) + coverage);
        }
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

    private void printGraph(DirectedGraph<String, Link> g, File f, String leftFlank, String rightFlank) {
        try {
            PrintStream o = new PrintStream(f);

            o.println("digraph G {");
            o.println("    graph [fontname = \"Courier\"];");
            o.println("    node [fontname = \"Courier\"];");
            o.println("    edge [fontname = \"Courier\"];");

            for (String v : g.vertexSet()) {
                if (v.equals(leftFlank) || v.equals(rightFlank)) {
                    o.println("    \"" + v + "\" [style=\"filled\" fillcolor=\"red\"];");
                } else {
                    o.println("    \"" + v + "\";");
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
                    o.println("    \"" + s + "\" -> \"" + t + "\" [ penwidth=" + Math.log10(cov + 1) + " color=\"" + col + "\" ];");
                }
            }

            o.println("}");

            o.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    private Set<String> walk(CortexGraph cg, CortexLinksMap linksMap, int color, String from, String to) {
        Set<String> reconstructionKmers = new HashSet<String>();
        StringBuilder reconstructionFwd = new StringBuilder();
        StringBuilder reconstructionRev = new StringBuilder();

        reconstructionKmers.add(from);
        reconstructionKmers.add(to);

        String thisKmer = from;
        reconstructionFwd.append(thisKmer);
        while (thisKmer != null) {
            Collection<String> kmers = CortexUtils.getPrevKmers(cg, thisKmer, color);

            if (kmers.size() == 1) {
                String tKmer = kmers.iterator().next();

                reconstructionKmers.add(tKmer);
                reconstructionFwd.insert(0, tKmer.charAt(0));

                thisKmer = tKmer;
            } else {
                for (String tKmer : kmers) { reconstructionKmers.add(tKmer); }
                thisKmer = null;
            }
        }

        thisKmer = from;
        while (thisKmer != null) {
            Collection<String> kmers = CortexUtils.getNextKmers(cg, thisKmer, color);

            if (kmers.size() == 1) {
                String tKmer = kmers.iterator().next();

                reconstructionKmers.add(tKmer);
                reconstructionFwd.append(tKmer.charAt(tKmer.length() - 1));

                thisKmer = tKmer;
            } else {
                for (String tKmer : kmers) { reconstructionKmers.add(tKmer); }
                thisKmer = null;
            }
        }

        thisKmer = to;
        reconstructionRev.append(thisKmer);
        while (thisKmer != null) {
            Collection<String> kmers = CortexUtils.getPrevKmers(cg, thisKmer, color);

            if (kmers.size() == 1) {
                String tKmer = kmers.iterator().next();

                reconstructionKmers.add(tKmer);
                reconstructionRev.insert(0, tKmer.charAt(0));

                thisKmer = tKmer;
            } else {
                for (String tKmer : kmers) { reconstructionKmers.add(tKmer); }
                thisKmer = null;
            }
        }

        thisKmer = to;
        while (thisKmer != null) {
            Collection<String> kmers = CortexUtils.getNextKmers(cg, thisKmer, color);

            if (kmers.size() == 1) {
                String tKmer = kmers.iterator().next();

                reconstructionKmers.add(tKmer);
                reconstructionRev.append(tKmer.charAt(tKmer.length() - 1));

                thisKmer = tKmer;
            } else {
                for (String tKmer : kmers) { reconstructionKmers.add(tKmer); }
                thisKmer = null;
            }
        }

        Set<String> alleles = new HashSet<String>();
        if (linksMap.getCortexLinks().hasColor(color)) {
            for (String reconstructionKmerFw : reconstructionKmers) {
                CortexKmer kmer = new CortexKmer(reconstructionKmerFw);
                if (linksMap.containsKey(kmer)) {
                    for (CortexJunctionsRecord cjr : linksMap.get(kmer).getJunctions()) {
                        String reconstruction = cjr.getSeq();

                        if (kmer.isFlipped() != cjr.isFlipped()) {
                            reconstruction = SequenceUtils.reverseComplement(reconstruction);
                        }

                        if (reconstruction.contains(from) && reconstruction.contains(to)) {
                            int fromIndex = reconstruction.indexOf(from) + from.length();
                            int toIndex = reconstruction.lastIndexOf(to);

                            String allele = reconstruction.substring(fromIndex, toIndex);
                            alleles.add(allele);
                        }
                    }
                }
            }
        } else {
            Set<String> reconstructions = new HashSet<String>();
            reconstructions.add(reconstructionFwd.toString());
            reconstructions.add(reconstructionRev.toString());

            for (String reconstruction : reconstructions) {
                if (reconstruction.contains(from) && reconstruction.contains(to)) {
                    int fromIndex = reconstruction.indexOf(from) + from.length();
                    int toIndex = reconstruction.lastIndexOf(to);

                    String allele = reconstruction.substring(fromIndex, toIndex);
                    alleles.add(allele);
                }
            }
        }

        return alleles;
    }

    @Override
    public void execute() {
        Set<CortexKmer> novelKmers = loadNovelKmers();
        int kmerSize = novelKmers.iterator().next().length();

        //CortexLinksMap linksMap = new CortexLinksMap(LINKS);
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
                printGraph(g, out, leftFlank, rightFlank);

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
                printGraph(g, out, leftFlank, rightFlank);

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
                printGraph(g, out, leftFlank, rightFlank);

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
                printGraph(g, out, leftFlank, rightFlank);

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
                    log.debug("     : {}",   childSequenceLeft);
                    log.debug("     : {}",   childSequenceRight);

                    Set<String> vertices = new HashSet<String>(g.vertexSet());

                    for (String sk : vertices) {
                        CortexKmer ck = new CortexKmer(sk);

                        if (linksMap.get(GRAPH.getColor(c).getSampleName()).containsKey(ck)) {
                            for (CortexJunctionsRecord cjr : linksMap.get(GRAPH.getColor(c).getSampleName()).get(ck).getJunctions()) {
                                String linkSeq = cjr.getSeq();

                                if (ck.isFlipped() != cjr.isFlipped()) {
                                    linkSeq = SequenceUtils.reverseComplement(linkSeq);
                                }

                                if (linkSeq.contains(leftFlank) || linkSeq.contains(rightFlank)) {
                                    for (int q = 0; q < linkSeq.length() - kmerSize; q++) {
                                        String tk = linkSeq.substring(q, q + kmerSize);
                                        String nk = linkSeq.substring(q + 1, q + 1 + kmerSize);
                                        g.addVertex(tk);
                                        g.addVertex(nk);

                                        if (!g.containsEdge(tk, nk)) {
                                            g.addEdge(tk, nk, new Link());
                                        } else {
                                            g.getEdge(tk, nk).addCoverage(c, cjr.getCoverage(0));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                log.info("{} {} {}:{}-{} {}", count, full, chr, start, stop, eventType);

                printGraph(g, out, leftFlank, rightFlank);
                break;
            }
        }

        log.info("  {}", Joiner.on("; ").withKeyValueSeparator("=").join(eventTypes));
    }
}
