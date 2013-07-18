package uk.ac.ox.well.indiana.utils.assembly;

import org.jgrapht.DirectedGraph;
import org.jgrapht.Graph;
import org.jgrapht.Graphs;
import org.jgrapht.experimental.dag.DirectedAcyclicGraph;
import org.jgrapht.ext.VertexNameProvider;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleDirectedGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexMap;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

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
        soughtKmers.add(new CortexKmer("TATATATATATATATATTATATGTGTATATG"));
        soughtKmers.add(new CortexKmer("CATATACACATATAATATATATATATATATA"));

        if (cortexMap.containsKey(panelKmer)) {
            if (soughtKmers.contains(panelKmer)) {
                System.out.println("Saw sought kmer '" + panelKmer + "'");
            }

            Set<CortexKmer> seenKmers = new HashSet<CortexKmer>();

            // First, add our panel kmer to the graph
            g.addVertex(panelKmer);

            // Next, figure out what the next edges are
            CortexRecord record = cortexMap.get(panelKmer);
            Collection<Byte> edges = goLeft ? record.getInEdgesAsBytes(color) : record.getOutEdgesAsBytes(color);
            byte[] currentKmer = record.getKmerAsBytes();
            //CortexKmer currentCortexKmer = panelKmer;

            boolean first = true;

            // If there's only one next edge, proceed
            while (record != null && !seenKmers.contains(record.getKmer()) && edges.size() == 1) {
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

                record = cortexMap.get(nextKmer);

                // Check to make sure this next kmer actually exists
                if (record != null) {
                    // Now add it to the graph
                    g.addVertex(nextKmer);
                    if (first) {
                        g.addEdge(panelKmer, nextKmer);
                        first = false;
                    } else {
                        g.addEdge(new CortexKmer(currentKmer), nextKmer);
                    }

                    if (goLeft) {
                        edges = nextKmer.isFlipped() ? record.getOutEdgesComplementAsBytes(color) : record.getInEdgesAsBytes(color);
                    } else {
                        edges = nextKmer.isFlipped() ? record.getInEdgesComplementAsBytes(color) : record.getOutEdgesAsBytes(color);
                    }
                    currentKmer = newKmer;
                    //currentCortexKmer = record.getKmer();
                }
            }

            // If we're allowing ourselves to navigate
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

                    g.addVertex(k2);
                    g.addEdge(ck, k2);

                    if (soughtKmers.contains(k2)) {
                        System.out.println("Saw sought kmer '" + k2 + "'");
                    }

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

        buildLocalGraph(color, pk, true, maxForksLeft, g);
        buildLocalGraph(color, pk, false, maxForksRight, g);

        return g;
    }

    public DirectedGraph<CortexKmer, DefaultEdge> buildLocalGraph(int color, CortexKmer panelKmer, int maxForks) {
        return buildLocalGraph(color, panelKmer, maxForks, maxForks);
    }
}
