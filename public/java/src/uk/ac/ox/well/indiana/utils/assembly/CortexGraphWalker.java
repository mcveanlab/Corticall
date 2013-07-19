package uk.ac.ox.well.indiana.utils.assembly;

import org.jgrapht.DirectedGraph;
import org.jgrapht.Graph;
import org.jgrapht.Graphs;
import org.jgrapht.experimental.dag.DirectedAcyclicGraph;
import org.jgrapht.ext.DOTExporter;
import org.jgrapht.ext.VertexNameProvider;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleDirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexMap;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class CortexGraphWalker {
    private CortexMap cortexMap;

    public CortexGraphWalker(CortexMap cortexMap) {
        this.cortexMap = cortexMap;
    }

    public CortexKmer buildContig(int color, CortexKmer panelKmer) {
        if (cortexMap.containsKey(panelKmer)) {
            StringBuilder contig = new StringBuilder();

            Set<CortexKmer> seenKmers = new HashSet<CortexKmer>();

            // start with left edges
            CortexRecord record = cortexMap.get(panelKmer);
            Collection<Byte> leftEdges = record.getInEdgesAsBytes(color);
            byte[] currentKmer = record.getKmerAsBytes();

            contig.append(new String(currentKmer));

            while (record != null && !seenKmers.contains(record.getKmer()) && leftEdges.size() == 1) {
                seenKmers.add(record.getKmer());

                byte[] newKmer = new byte[currentKmer.length];

                newKmer[0] = leftEdges.iterator().next();
                System.arraycopy(currentKmer, 0, newKmer, 1, currentKmer.length - 1);

                contig.insert(0, SequenceUtils.nucleotideByteToString(newKmer[0]));

                CortexKmer leftKmer = new CortexKmer(newKmer);
                record = cortexMap.get(leftKmer);

                if (record != null) {
                    leftEdges = leftKmer.isFlipped() ? record.getOutEdgesComplementAsBytes(color) : record.getInEdgesAsBytes(color);
                    currentKmer = newKmer;
                }
            }

            // now do right edges
            record = cortexMap.get(panelKmer);
            Collection<Byte> rightEdges = record.getOutEdgesAsBytes(color);
            currentKmer = record.getKmerAsBytes();

            seenKmers.remove(panelKmer);

            while (record != null && !seenKmers.contains(record.getKmer()) && rightEdges.size() == 1) {
                seenKmers.add(record.getKmer());

                byte[] newKmer = new byte[currentKmer.length];

                System.arraycopy(currentKmer, 1, newKmer, 0, currentKmer.length - 1);
                newKmer[currentKmer.length - 1] = rightEdges.iterator().next();

                contig.append(SequenceUtils.nucleotideByteToString(newKmer[currentKmer.length - 1]));

                CortexKmer rightKmer = new CortexKmer(newKmer);
                record = cortexMap.get(rightKmer);

                if (record != null) {
                    rightEdges = rightKmer.isFlipped() ? record.getInEdgesComplementAsBytes(color) : record.getOutEdgesAsBytes(color);
                    currentKmer = newKmer;
                }
            }

            return new CortexKmer(contig.toString());
        }

        return null;
    }

    public Map<CortexKmer, Set<CortexKmer>> buildContigs(int color, Set<CortexKmer> panelKmers) {
        Set<CortexKmer> usedKmers = new HashSet<CortexKmer>();

        Map<CortexKmer, Set<CortexKmer>> contigs = new HashMap<CortexKmer, Set<CortexKmer>>();

        for (CortexKmer panelKmer : panelKmers) {
            if (!usedKmers.contains(panelKmer)) {
                CortexKmer contig = buildContig(color, panelKmer);

                if (contig != null) {
                    Set<CortexKmer> seedKmers = new HashSet<CortexKmer>();

                    byte[] contigBytes = contig.getKmerAsBytes();

                    for (int i = 0; i <= contigBytes.length - panelKmer.length(); i++) {
                        byte[] kmer = new byte[panelKmer.length()];
                        System.arraycopy(contigBytes, i, kmer, 0, panelKmer.length());

                        CortexKmer ckmer = new CortexKmer(kmer);
                        if (panelKmers.contains(ckmer)) {
                            seedKmers.add(ckmer);
                        }
                    }

                    usedKmers.addAll(seedKmers);

                    contigs.put(contig, seedKmers);
                }
            }
        }

        return contigs;
    }

    private DirectedGraph<CortexKmer, DefaultEdge> buildLocalGraph(int color, CortexKmer panelKmer, boolean goLeft, int maxForks, DirectedGraph<CortexKmer, DefaultEdge> g) {
        Set<CortexKmer> soughtKmers = new HashSet<CortexKmer>();
        soughtKmers.add(new CortexKmer("TTTTTATTTAATAAATTTGTTTTTATTTTAT"));

        if (cortexMap.containsKey(panelKmer)) {
            System.out.println("Start: color=" + color + " kmer=" + panelKmer + " goLeft=" + goLeft + " maxForks=" + maxForks + " numVertices=" + g.vertexSet().size());

            if (soughtKmers.contains(panelKmer) && g.containsVertex(panelKmer)) {
                System.out.println("Found kmer '" + panelKmer + "'");
            }

            Set<CortexKmer> seenKmers = new HashSet<CortexKmer>();

            // First, add our panel kmer to the graph
            g.addVertex(panelKmer);

            // Next, figure out what the next edges are
            CortexRecord record = cortexMap.get(panelKmer);
            Collection<Byte> edges = goLeft ? record.getInEdgesAsBytes(color) : record.getOutEdgesAsBytes(color);
            byte[] currentKmer = record.getKmerAsBytes();

            boolean first = true;

            // If there's only one next edge, proceed
            while (record != null && !seenKmers.contains(record.getKmer()) && edges.size() == 1) {
                System.out.println("\tTraverse: color=" + color + " kmer=" + new CortexKmer(currentKmer) + " goLeft=" + goLeft + " maxForks=" + maxForks + " numVertices=" + g.vertexSet().size());

                seenKmers.add(record.getKmer());

                // Construct the next kmer based on this edge
                byte[] newKmer = new byte[currentKmer.length];

                if (goLeft) {
                    newKmer[0] = edges.iterator().next();
                    System.arraycopy(currentKmer, 0, newKmer, 1, currentKmer.length - 1);
                } else {
                    System.arraycopy(currentKmer, 1, newKmer, 0, currentKmer.length - 1);
                    newKmer[currentKmer.length - 1] = edges.iterator().next();
                }

                CortexKmer nextKmer = new CortexKmer(newKmer);

                System.out.println("\tNextkmer: color=" + color + " kmer=" + nextKmer + " goLeft=" + goLeft + " maxForks=" + maxForks + " numVertices=" + g.vertexSet().size());

                if (soughtKmers.contains(nextKmer) && g.containsVertex(nextKmer)) {
                    System.out.println("Found kmer '" + nextKmer + "'");
                }

                record = cortexMap.get(nextKmer);

                // Check to make sure this next kmer actually exists
                if (record != null) {
                    // Now add it to the graph
                    g.addVertex(nextKmer);
                    if (first) {
                        if (goLeft) {
                            g.addEdge(panelKmer, nextKmer);
                        } else {
                            g.addEdge(nextKmer, panelKmer);
                        }
                        first = false;
                    } else {
                        if (goLeft) {
                            g.addEdge(new CortexKmer(currentKmer), nextKmer);
                        } else {
                            g.addEdge(nextKmer, new CortexKmer(currentKmer));
                        }
                    }

                    if (goLeft) {
                        edges = nextKmer.isFlipped() ? record.getOutEdgesComplementAsBytes(color) : record.getInEdgesAsBytes(color);
                    } else {
                        edges = nextKmer.isFlipped() ? record.getInEdgesComplementAsBytes(color) : record.getOutEdgesAsBytes(color);
                    }
                    currentKmer = newKmer;
                }

                /*
                DOTExporter<CortexKmer, DefaultEdge> exporter = new DOTExporter<CortexKmer, DefaultEdge>(new CortexKmerNameProvider(), new CortexKmerNameProvider(), null);
                try {
                    exporter.export(new FileWriter("test9.dot"), g);
                    Thread.sleep(100);
                } catch (IOException e) {
                    e.printStackTrace();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
                */
            }

            // If we're allowing ourselves to navigate past forks, handle it here
            if (maxForks > 0 && edges.size() > 1) {
                CortexKmer ck = new CortexKmer(currentKmer);

                for (Byte edge : edges) {
                    byte[] newKmer = new byte[currentKmer.length];

                    if (goLeft) {
                        newKmer[0] = edge;
                        System.arraycopy(currentKmer, 0, newKmer, 1, currentKmer.length - 1);
                    } else {
                        System.arraycopy(currentKmer, 1, newKmer, 0, currentKmer.length - 1);
                        newKmer[currentKmer.length - 1] = edge;
                    }

                    CortexKmer k2 = new CortexKmer(newKmer);

                    if (soughtKmers.contains(k2) && g.containsVertex(k2)) {
                        System.out.println("Found kmer '" + k2 + "'");
                    }

                    g.addVertex(k2);
                    if (goLeft) {
                        g.addEdge(ck, k2);
                    } else {
                        g.addEdge(k2, ck);
                    }

                    System.out.println("Branch from: color=" + color + " kmer=" + ck + " goLeft=" + goLeft + " maxForks=" + maxForks + " numVertices=" + g.vertexSet().size());
                    System.out.println("Branch to  : color=" + color + " kmer=" + k2 + " goLeft=" + goLeft + " maxForks=" + maxForks + " numVertices=" + g.vertexSet().size());

                    g = buildLocalGraph(color, k2, k2.isFlipped() ? !goLeft : goLeft, maxForks - 1, g);
                }
            }

            return g;
        }

        return null;
    }

    public DirectedGraph<CortexKmer, DefaultEdge> buildLocalGraph(int color, CortexKmer panelKmer, int maxForksLeft, int maxForksRight) {
        DirectedGraph<CortexKmer, DefaultEdge> g = new DefaultDirectedGraph<CortexKmer, DefaultEdge>(DefaultEdge.class);

        // This is here to address a weird bug.  If you call this method with a CortexKmer that was constructed with
        // non-alphanumerically-lowest kmer, the graph object will contain one extra vertex with the non-alphanumerically
        // lowest kmer.  This is solved by simply ensuring the CortexKmer that we start with does *not* have the
        // isFlipped == true attribute.
        CortexKmer pk = new CortexKmer(panelKmer.getKmerAsString());

        g = buildLocalGraph(color, pk, true, maxForksLeft, g);
        g = buildLocalGraph(color, pk, false, maxForksRight, g);

        return g;
    }

    public DirectedGraph<CortexKmer, DefaultEdge> buildLocalGraph(int color, CortexKmer panelKmer, int maxForks) {
        return buildLocalGraph(color, panelKmer, maxForks, maxForks);
    }

    public DirectedGraph<CortexKmer, DefaultEdge> buildLocalGraph(CortexKmer panelKmer, int maxForksLeft, int maxForksRight) {
        DirectedGraph<CortexKmer, DefaultEdge> g = new DefaultDirectedGraph<CortexKmer, DefaultEdge>(DefaultEdge.class);
        CortexKmer pk = new CortexKmer(panelKmer.getKmerAsString());

        for (int color = 0; color < cortexMap.getGraph().getNumColors(); color++) {
            g = buildLocalGraph(color, pk, true, maxForksLeft, g);
            g = buildLocalGraph(color, pk, false, maxForksRight, g);
        }

        return g;
    }

    public DirectedGraph<CortexKmer, DefaultEdge> buildLocalGraph(CortexKmer panelKmer, int maxForks) {
        return buildLocalGraph(panelKmer, maxForks, maxForks);
    }
}
