package uk.ac.ox.well.indiana.commands.prg;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

@Description(text="Create a population reference graph from the given FASTA file")
public class createprg extends Module {
    @Argument(fullName="fasta", shortName="f", doc="FASTA file")
    public FastaSequenceFile FASTA;

    @Argument(fullName="proteinKmerSize", shortName="pk", doc="Kmer size (in amino-acid space)")
    public Integer PROTEIN_KMER_SIZE = 7;

    @Output
    public File out;

    private void writeGraph(DirectedGraph<String, DefaultEdge> g, File o) {
        try {
            PrintStream ps = new PrintStream(o);

            String indent = "  ";

            ps.println("digraph G {");
            ps.println(indent + "rankdir=\"LR\";");
            ps.println(indent + "edge [ dir=both arrowhead=none arrowtail=none ];");
            ps.println(indent + "node [ shape=none fontname=courier fontsize=9 ];");

            for (String vertex : g.vertexSet()) {
                //ps.println(indent + vertex + " [ shape=none label=\"\" ];");
                ps.println(indent + vertex);
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
        // Construct initial de Bruijn graph
        log.info("Constructing initial de Bruijn graph...");
        int kmerSize = 3*PROTEIN_KMER_SIZE;

        DirectedGraph<String, DefaultEdge> dbg = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);

        Map<String, String> sequences = new HashMap<String, String>();
        ReferenceSequence rseq;
        while ((rseq = FASTA.nextSequence()) != null) {
            String pkmer = null;
            if (rseq.getName().endsWith(".gene")) {
                String seq = new String(rseq.getBases());

                for (int i = 0; i <= seq.length() - kmerSize; i++) {
                    String ckmer = seq.substring(i, i + kmerSize);
                    ckmer = ckmer.replace("*", "_");

                    dbg.addVertex(ckmer);

                    if (pkmer != null) {
                        dbg.addEdge(pkmer, ckmer);
                    }

                    pkmer = ckmer;
                }

                sequences.put(rseq.getName(), seq);
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

        for (String seqName : sequences.keySet()) {
            String seq = sequences.get(seqName);

            StringBuilder contig = new StringBuilder();
            String prevContig = null;

            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                String kmer = seq.substring(i, i + kmerSize);

                if (contig.length() == 0) {
                    contig.append(kmer);
                } else {
                    contig.append(kmer.charAt(kmer.length() - 1));
                }

                if (branchVertices.contains(kmer) && i > 0) {
                    sg.addVertex(contig.toString());

                    if (prevContig != null) {
                        sg.addEdge(prevContig, contig.toString());
                    }

                    prevContig = contig.toString();
                    contig = new StringBuilder();
                    contig.append(kmer);
                }
            }
        }

        writeGraph(sg, out);
    }
}
