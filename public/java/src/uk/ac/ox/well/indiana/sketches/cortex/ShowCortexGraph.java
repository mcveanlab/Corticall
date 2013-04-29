package uk.ac.ox.well.indiana.sketches.cortex;

import org.jgrapht.DirectedGraph;
import org.jgrapht.ext.ComponentAttributeProvider;
import org.jgrapht.ext.DOTExporter;
import org.jgrapht.ext.StringNameProvider;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class ShowCortexGraph extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    //@Argument(fullName="colors", shortName="c", doc="Colors to process")
    //public ArrayList<Integer> COLORS;

    @Output
    //public PrintStream out;
    public String out;

    private class CortexAttributeProvider implements ComponentAttributeProvider<String> {
        private HashMap<String, Map<String, String>> attributes = new HashMap<String, Map<String, String>>();

        public void setComponentAttributes(String s, String key, String value) {
            if (!attributes.containsKey(s)) {
                attributes.put(s, new HashMap<String, String>());
            }

            attributes.get(s).put(key, value);
        }

        @Override
        public Map<String, String> getComponentAttributes(String s) {
            return attributes.get(s);
        }
    }

    private void simplifyGraph(DirectedGraph<String, DefaultEdge> directedGraph, String kmer) {
        HashSet<DefaultEdge> seenEdges = new HashSet<DefaultEdge>();

        simplifyGraph(directedGraph, kmer, seenEdges, 0);
    }

    private void simplifyGraph(DirectedGraph<String, DefaultEdge> directedGraph, String kmer, HashSet<DefaultEdge> seenEdges, int indentLevel) {
        String indent = "";
        for (int i = 0; i < indentLevel; i++) {
            indent += " ";
        }

        log.info(indent + kmer);

        Set<DefaultEdge> edges = directedGraph.outgoingEdgesOf(kmer);

        for (DefaultEdge edge : edges) {
            if (!seenEdges.contains(edge)) {
                //seenEdges.add(edge);

                String target = directedGraph.getEdgeTarget(edge);

                int indentAdd = edges.size() == 1 ? 0 : 1;

                simplifyGraph(directedGraph, target, seenEdges, indentLevel + indentAdd);
            }
        }
    }

    public void execute() {
        DirectedGraph<String, DefaultEdge> directedGraph = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
        CortexAttributeProvider cap = new CortexAttributeProvider();

        //String startKmer = null;

        int count = 0;

        for (CortexRecord cr : CORTEX_GRAPH) {
            count++;
            if (count % 1000000 == 0) {
                log.info("Processed {}/{} records", count, CORTEX_GRAPH.getNumRecords());
            }

            String kmer = cr.getKmerAsString();

            directedGraph.addVertex(kmer);

            //boolean isStartKmer = true;

            for (int color = 0; color < CORTEX_GRAPH.getNumColors(); color++) {
                String edges = (cr.getEdgeAsStrings()[color]).toUpperCase();

                for (int inEdgeIndex = 0; inEdgeIndex < 4; inEdgeIndex++) {
                    if (edges.charAt(inEdgeIndex) != '.') {
                        String inEdgeKmer = SequenceUtils.alphanumericallyLowestOrientation(edges.charAt(inEdgeIndex) + kmer.substring(0, kmer.length() - 1));

                        if (!directedGraph.containsVertex(inEdgeKmer)) {
                            directedGraph.addVertex(inEdgeKmer);
                        }

                        directedGraph.addEdge(inEdgeKmer, kmer);

                        //isStartKmer = false;
                    }
                }

                for (int outEdgeIndex = 4; outEdgeIndex < 8; outEdgeIndex++) {
                    if (edges.charAt(outEdgeIndex) != '.') {
                        String outEdgeKmer = SequenceUtils.alphanumericallyLowestOrientation(kmer.substring(1, kmer.length()) + edges.charAt(outEdgeIndex));

                        if (!directedGraph.containsVertex(outEdgeKmer)) {
                            directedGraph.addVertex(outEdgeKmer);
                        }

                        directedGraph.addEdge(kmer, outEdgeKmer);
                    }
                }
            }

            int[] coverages = cr.getCoverages();
            if (coverages[0] > 0 && coverages[1] > 0) {
                cap.setComponentAttributes(kmer, "label", "both");
            } else if (coverages[0] > 0 && coverages[1] == 0) {
                cap.setComponentAttributes(kmer, "label", "s1");
            } else {
                cap.setComponentAttributes(kmer, "label", "s2");
            }

            //if (isStartKmer) {
                //log.info("start kmer: {}", cr);

                //startKmer = kmer;
            //}
        }

        /*
        DirectedGraph<String, DefaultEdge> dg = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);

        AACTGGTCA
        AACTCGTCA

        AAC
         ACT
          CTG
           TGG
            GGT
             GTC
              TCA

         AAC
          ACT
           CTC
            TCG
             CGT
              GTC
               TCA

        dg.addVertex("AAC");
        dg.addVertex("ACT");
        dg.addVertex("CTG");
        dg.addVertex("TGG");
        dg.addVertex("GGT");
        dg.addVertex("GTC");
        dg.addVertex("TCA");

        dg.addVertex("CTC");
        dg.addVertex("TCG");
        dg.addVertex("CGT");

        dg.addEdge("AAC", "ACT");
        dg.addEdge("ACT", "CTG");
        dg.addEdge("CTG", "TGG");
        dg.addEdge("TGG", "GGT");
        dg.addEdge("GGT", "GTC");
        dg.addEdge("GTC", "TCA");

        dg.addEdge("ACT", "CTC");
        dg.addEdge("CTC", "TCG");
        dg.addEdge("TCG", "CGT");
        dg.addEdge("CGT", "GTC");

        simplifyGraph(dg, "AAC");
        */

        /*
        if (startKmer != null) {
            simplifyGraph(directedGraph, startKmer);
        }
        */

        DOTExporter<String, DefaultEdge> exporter = new DOTExporter<String, DefaultEdge>(new StringNameProvider<String>(), null, null, cap, null);

        try {
            //exporter.export(new FileWriter(out), dg);
            exporter.export(new FileWriter(out), directedGraph);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
