package uk.ac.ox.well.cortexjdk.utils.assembler;

import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

/**
 * Created by kiran on 15/09/2017.
 */
public class TempLinkThreader {
    private TempLinkThreader() {}

    public static CortexLinks buildLinks(CortexGraph graph, Map<String, Collection<String>> haplotypeLists, String sampleName) {
        int color = graph.getColorForSampleName(sampleName);

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

        Collection<String> haplotypes = haplotypeLists.get(sampleName);


        File tempFile = null;
        try {
            tempFile = File.createTempFile("tempgraph", ".ctp.bgz");
        } catch (IOException e) {
            throw new CortexJDKException("Could not get a temp file for links creation");
        }
        tempFile.deleteOnExit();



        return null;
    }
}
