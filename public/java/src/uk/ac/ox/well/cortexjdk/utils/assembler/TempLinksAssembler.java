package uk.ac.ox.well.cortexjdk.utils.assembler;

import com.google.common.base.Joiner;
import org.jetbrains.annotations.NotNull;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.util.*;

/**
 * Created by kiran on 15/09/2017.
 */
public class TempLinksAssembler {
    private TempLinksAssembler() {}

    public static CortexLinks buildLinks(CortexGraph graph, Map<String, Collection<String>> haplotypeLists, String sampleName) {
        DirectedGraph<String, DefaultEdge> g = loadGraph(graph, graph.getColorForSampleName(sampleName));

        for (String haplotypeFwd : haplotypeLists.get(sampleName)) {
            String haplotypeRev = SequenceUtils.reverseComplement(haplotypeFwd);

            for (String haplotype : Arrays.asList(haplotypeFwd, haplotypeRev)) {
                List<List<String>> links = new ArrayList<>();
                List<List<String>> pos = new ArrayList<>();

                for (int i = 0; i <= haplotype.length() - graph.getKmerSize() - 1; i++) {
                    String sk0 = haplotype.substring(i, i + graph.getKmerSize());
                    String sk1 = haplotype.substring(i + 1, i + graph.getKmerSize() + 1);
                    String edge = haplotype.substring(i + graph.getKmerSize(), i + graph.getKmerSize() + 1);

                    if (g.outDegreeOf(sk0) > 1 && g.containsVertex(sk1)) {
                        links.add(new ArrayList<>());
                        pos.add(new ArrayList<>());

                        for (int l = 0; l < links.size(); l++) {
                            links.get(l).add(edge);

                            for (int j = i - 1; j >= 0; j--) {
                                String skj = haplotype.substring(j, j + graph.getKmerSize());
                                if (g.inDegreeOf(skj) <= 1 && g.outDegreeOf(skj) <= 1) {
                                    pos.get(l).add(skj);
                                }
                            }
                        }
                    }
                }

                System.out.println(haplotype);
                for (int l = 0; l < links.size(); l++) {
                    System.out.println("  " + l + " " + pos.get(l).get(0) + " " + Joiner.on("").join(links.get(l)));
                }
            }
        }

        /*
        File tempFile = null;
        try {
            tempFile = File.createTempFile("tempgraph", ".ctp.bgz");
        } catch (IOException e) {
            throw new CortexJDKException("Could not get a temp file for links creation");
        }
        tempFile.deleteOnExit();
        */

        return null;
    }

    @NotNull
    private static DirectedGraph<String, DefaultEdge> loadGraph(CortexGraph graph, int color) {
        DirectedGraph<String, DefaultEdge> g = new DefaultDirectedGraph<>(DefaultEdge.class);

        for (CortexRecord cr : graph) {
            if (cr.getCoverage(color) > 0) {
                String fwd = cr.getKmerAsString();
                g.addVertex(fwd);

                Collection<String> inEdges = cr.getInEdgesAsStrings(color, false);
                for (String inEdge : inEdges) {
                    String inSk = inEdge + fwd.substring(0, fwd.length() - 1);
                    g.addVertex(inSk);
                    g.addEdge(inSk, fwd);
                }

                Collection<String> outEdges = cr.getOutEdgesAsStrings(color, false);
                for (String outEdge : outEdges) {
                    String outSk = fwd.substring(1, fwd.length()) + outEdge;
                    g.addVertex(outSk);
                    g.addEdge(fwd, outSk);
                }

                String rev = SequenceUtils.reverseComplement(cr.getKmerAsString());
                g.addVertex(rev);

                inEdges = cr.getOutEdgesAsStrings(color, true);
                for (String inEdge : inEdges) {
                    String inSk = inEdge + rev.substring(0, rev.length() - 1);
                    g.addVertex(inSk);
                    g.addEdge(inSk, rev);
                }

                outEdges = cr.getInEdgesAsStrings(color, true);
                for (String outEdge : outEdges) {
                    String outSk = rev.substring(1, rev.length()) + outEdge;
                    g.addVertex(outSk);
                    g.addEdge(rev, outSk);
                }
            }
        }
        return g;
    }
}
