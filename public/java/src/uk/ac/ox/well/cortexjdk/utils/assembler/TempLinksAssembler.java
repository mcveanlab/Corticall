package uk.ac.ox.well.cortexjdk.utils.assembler;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import org.jetbrains.annotations.NotNull;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.json.JSONArray;
import org.json.JSONObject;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexJunctionsRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinksRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.*;
import java.util.*;

/**
 * Created by kiran on 15/09/2017.
 */
public class TempLinksAssembler {
    private TempLinksAssembler() {}

    public static CortexLinks buildLinks(CortexGraph graph, Map<String, Collection<String>> haplotypeLists, String sampleName) {
        DirectedGraph<String, DefaultEdge> g = loadGraph(graph, graph.getColorForSampleName(sampleName));

        Map<CanonicalKmer, Set<CortexJunctionsRecord>> linkMap = new HashMap<>();

        for (String haplotypeFwd : haplotypeLists.get(sampleName)) {
            String haplotypeRev = SequenceUtils.reverseComplement(haplotypeFwd);

            for (String haplotype : Arrays.asList(haplotypeFwd, haplotypeRev)) {
                Map<String, String> links = new HashMap<>();

                for (int j = 0; j <= haplotype.length() - graph.getKmerSize() - 1; j++) {
                    String sk0 = haplotype.substring(j, j + graph.getKmerSize());
                    String sk1 = haplotype.substring(j + 1, j + graph.getKmerSize() + 1);
                    String edge = haplotype.substring(j + graph.getKmerSize(), j + graph.getKmerSize() + 1);

                    if (g.outDegreeOf(sk0) > 1 && g.containsVertex(sk1)) {
                        for (int i = 1; i <= j; i++) {
                            String ski = haplotype.substring(i, i + graph.getKmerSize());
                            if (g.inDegreeOf(ski) > 1) {
                                String skim1 = haplotype.substring(i-1, i-1 + graph.getKmerSize());

                                if (!links.containsKey(skim1)) {
                                    links.put(skim1, "");
                                }

                                links.put(skim1, links.get(skim1) + edge);
                            }
                        }
                    }
                }

                for (String anchor : links.keySet()) {
                    CanonicalKmer ck = new CanonicalKmer(anchor);

                    if (!linkMap.containsKey(ck)) {
                        linkMap.put(ck, new HashSet<>());
                    }

                    linkMap.get(ck).add(new CortexJunctionsRecord(!ck.isFlipped(), links.get(anchor).length(), links.get(anchor).length(), new int[] { 1 }, links.get(anchor)));
                }
            }
        }

        int numPaths = 0;
        for (CanonicalKmer ck : linkMap.keySet()) {
            numPaths += linkMap.get(ck).size();
        }

        File tempFile;
        try {
            tempFile = File.createTempFile("templinks", ".ctp.bgz");
            tempFile.deleteOnExit();

            BlockCompressedOutputStream os = new BlockCompressedOutputStream(tempFile);

            os.write(constructLinksHeader(graph.getKmerSize(), graph.getNumRecords(), sampleName, linkMap.keySet().size(), numPaths).toString(8).getBytes());

            os.write("\n\n".getBytes());

            for (CanonicalKmer ck : linkMap.keySet()) {
                CortexLinksRecord clr = new CortexLinksRecord(ck.getKmerAsString(), new ArrayList<>(linkMap.get(ck)));

                os.write(clr.toString().getBytes());
                os.write("\n".getBytes());
            }

            os.write("\n".getBytes());

            os.close();

            return new CortexLinks(tempFile);
        } catch (IOException e) {
            throw new CortexJDKException("Could not get a temp file for links creation");
        }
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

    private static JSONObject constructLinksHeader(int kmerSize, long numKmersInGraph, String sample, int numKmersWithLinks, int numPaths) {
        JSONObject header = new JSONObject();
        header.put("file_format", "ctp");
        header.put("format_version", 4);
        header.put("file_key", 0);

        JSONArray colors = new JSONArray();
        JSONObject color = new JSONObject();
        color.put("colour", 0);
        color.put("sample", sample);
        color.put("total_sequence", 0);
        color.put("cleaned_tips", false);
        color.put("cleaned_unitigs", false);
        colors.put(color);

        JSONObject jsonGraph = new JSONObject();
        jsonGraph.put("num_colours", 1);
        jsonGraph.put("kmer_size", kmerSize);
        jsonGraph.put("num_kmers_in_graph", numKmersInGraph);
        jsonGraph.put("colours", colors);

        header.put("graph", jsonGraph);

        JSONObject jsonPaths = new JSONObject();
        jsonPaths.put("num_kmers_with_paths", numKmersWithLinks);
        jsonPaths.put("num_paths", numPaths);
        jsonPaths.put("path_bytes", numPaths);

        header.put("paths", jsonPaths);

        return header;
    }
}
