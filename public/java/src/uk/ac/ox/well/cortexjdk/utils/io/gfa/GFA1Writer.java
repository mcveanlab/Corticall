package uk.ac.ox.well.cortexjdk.utils.io.gfa;

import com.google.common.base.Joiner;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

public class GFA1Writer {
    public static void write(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, File dir, String name) {
        if (g.vertexSet().size() > 0) {
            try {
                PrintStream out = new PrintStream(dir.getAbsolutePath() + "/" + name + ".gfa");

                int kmerSize = g.vertexSet().iterator().next().getKmerAsString().length();

                out.println(Joiner.on("\t").join("H", "VN:Z:1.0"));

                int index = 0;
                Map<CortexVertex, Integer> indices = new HashMap<>();
                for (CortexVertex v : g.vertexSet()) {
                    out.println(Joiner.on("\t").join("S", index, v.getKmerAsString())); // , String.format("RC:i:%d", v.), String.format("AC:i:%d", avCov)));

                    indices.put(v, index);
                    index++;
                }

                for (CortexEdge e : g.edgeSet()) {
                    CortexVertex vs = g.getEdgeSource(e);
                    int vsid = indices.get(vs);
                    String vsStrand = vs.getKmerAsString().equalsIgnoreCase(vs.getCanonicalKmer().getKmerAsString()) ? "+" : "-";

                    CortexVertex vt = g.getEdgeTarget(e);
                    int vtid = indices.get(vt);
                    String vtStrand = vt.getKmerAsString().equalsIgnoreCase(vt.getCanonicalKmer().getKmerAsString()) ? "+" : "-";

                    out.println(Joiner.on("\t").join("L", vsid, vsStrand, vtid, vtStrand, String.format("%dM", kmerSize - 1)));
                }
            } catch (FileNotFoundException e) {
                throw new CortexJDKException("File not found", e);
            }
        }
    }
}
