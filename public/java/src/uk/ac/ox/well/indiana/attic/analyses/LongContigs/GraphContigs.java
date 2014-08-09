package uk.ac.ox.well.indiana.attic.analyses.LongContigs;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public class GraphContigs extends Module {
    @Argument(fullName="classifiedContigs", shortName="cc", doc="Classified contigs")
    public File CLASSIFIED_CONTIGS;

    @Argument(fullName="panelFasta", shortName="pf", doc="Panel")
    public FastaSequenceFile PANEL_FASTA;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Argument(fullName="select1", shortName="s1", doc="Select first ID to color")
    public String SELECT1;

    @Argument(fullName="select2", shortName="s2", doc="Select second ID to color")
    public String SELECT2;

    @Output
    public PrintStream out;

    private void writeGraph(DirectedGraph<String, DefaultEdge> g, File o, Map<String, Set<String>> ids) {
        try {
            PrintStream ps = new PrintStream(o);

            String indent = "  ";

            ps.println("digraph G {");
            ps.println(indent + "rankdir=\"LR\";");
            ps.println(indent + "edge [ dir=both arrowhead=none arrowtail=none ];");
            ps.println(indent + "node [ shape=none fontname=courier fontsize=9 ];");

            for (String vertex : g.vertexSet()) {
                String id = ids.containsKey(vertex) ? Joiner.on(",").join(ids.get(vertex)) : "";
                ps.println(indent + vertex + " [ label=\"" + id + "\" ];");
                //ps.println(indent + vertex);
            }

            for (DefaultEdge e : g.edgeSet()) {
                ps.println(indent + g.getEdgeSource(e) + " -> " + g.getEdgeTarget(e));
            }

            ps.println("}");
        } catch (FileNotFoundException e) {
            throw new IndianaException("Unable to write graph to file '" + o.getAbsolutePath() + "'", e);
        }
    }

    @Override
    public void execute() {
        Map<String, ReferenceSequence> panelSeqs = new HashMap<String, ReferenceSequence>();
        Map<CortexKmer, String> kmerToSeqMap = new HashMap<CortexKmer, String>();

        ReferenceSequence rseq;
        while ((rseq = PANEL_FASTA.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            panelSeqs.put(rseq.getName(), rseq);

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i , i + KMER_SIZE));

                kmerToSeqMap.put(kmer, rseq.getName());
            }
        }

        TableReader tr = new TableReader(CLASSIFIED_CONTIGS);

        for (Map<String, String> te : tr) {
            String name = te.get("name");
            String contig = te.get("contig");
            Integer kmersHB3 = Integer.valueOf(te.get(SELECT1));
            Integer kmers3D7 = Integer.valueOf(te.get(SELECT2));

            if (kmersHB3 > 0 && kmers3D7 > 0) {
                log.info("{} {} {}", name, kmersHB3, kmers3D7);

                Set<ReferenceSequence> seqs = new HashSet<ReferenceSequence>();
                Set<String> seqnames = new TreeSet<String>();

                for (int i = 0; i <= contig.length() - KMER_SIZE; i++) {
                    CortexKmer kmer = new CortexKmer(contig.substring(i, i + KMER_SIZE));

                    if (kmerToSeqMap.containsKey(kmer)) {
                        ReferenceSequence sseq = panelSeqs.get(kmerToSeqMap.get(kmer));
                        seqs.add(sseq);
                        seqnames.add(sseq.getName());
                    }
                }

                out.println(name + " " + Joiner.on(" ").join(seqnames));

                DirectedGraph<String, DefaultEdge> dbg = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
                Map<String, String> sequences = new HashMap<String, String>();

                sequences.put(name, contig);

                String pkmer = null;
                for (int i = 0; i <= contig.length() - KMER_SIZE; i++) {
                    String ckmer = contig.substring(i, i + KMER_SIZE);

                    dbg.addVertex(ckmer);

                    if (pkmer != null) {
                        dbg.addEdge(pkmer, ckmer);
                    }

                    pkmer = ckmer;
                }

                for (ReferenceSequence aseq : seqs) {
                    String seq = new String(aseq.getBases());

                    int fw = 0, rc = 0;
                    for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                        String kmerFw = seq.substring(i, i + KMER_SIZE);
                        String kmerRc = SequenceUtils.reverseComplement(kmerFw);

                        if      (dbg.containsVertex(kmerFw)) { fw++; }
                        else if (dbg.containsVertex(kmerRc)) { rc++; }
                    }

                    if (rc > fw) {
                        seq = SequenceUtils.reverseComplement(seq);
                    }

                    sequences.put(aseq.getName(), seq);

                    pkmer = null;
                    for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                        String ckmer = seq.substring(i, i + KMER_SIZE);

                        dbg.addVertex(ckmer);

                        if (pkmer != null) {
                            dbg.addEdge(pkmer, ckmer);
                        }

                        pkmer = ckmer;
                    }
                }

                // Find all of the branch points in the graph.
                log.info("Find branch kmers within the graph...");
                Set<String> branchVertices = new HashSet<String>();

                for (String vertex : dbg.vertexSet()) {
                    if (dbg.inDegreeOf(vertex) != 1 || dbg.outDegreeOf(vertex) != 1) {
                        branchVertices.add(vertex);
                    }
                }

                log.info("  Found {} branch kmers (kmers with more than one in or out degree)", branchVertices.size());

                // Rebuild the graph as a string graph
                log.info("Rebuilding the graph as a string graph...");

                DirectedGraph<String, DefaultEdge> sg = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
                Map<String, Set<String>> ids = new HashMap<String, Set<String>>();

                for (String seqName : sequences.keySet()) {
                    String seq = sequences.get(seqName);

                    StringBuilder acontig = new StringBuilder();
                    String prevContig = null;

                    for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                        String kmer = seq.substring(i, i + KMER_SIZE);

                        if (contig.length() == 0) {
                            acontig.append(kmer);
                        } else {
                            acontig.append(kmer.charAt(kmer.length() - 1));
                        }

                        if (branchVertices.contains(kmer) && i > 0) {
                            sg.addVertex(acontig.toString());

                            if (!ids.containsKey(acontig.toString())) {
                                ids.put(acontig.toString(), new TreeSet<String>());
                            }

                            ids.get(acontig.toString()).add(seqName);

                            if (prevContig != null) {
                                sg.addEdge(prevContig, acontig.toString());
                            }

                            prevContig = acontig.toString();
                            acontig = new StringBuilder();
                            acontig.append(kmer);
                        }
                    }
                }

                File o = new File(name + ".dot");
                writeGraph(sg, o, ids);

                try {
                    //Runtime.getRuntime().exec("dot -Tps2 " + o.getAbsolutePath() + " -o " + name + ".ps && ps2pdf " + name + ".ps");
                    Runtime.getRuntime().exec("dot -Tpng " + o.getAbsolutePath() + " -o " + name + ".png");
                } catch (IOException e) {
                    throw new IndianaException("Unable to convert .dot file to .pdf", e);
                }
            }
        }
    }
}
