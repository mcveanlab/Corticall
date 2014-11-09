package uk.ac.ox.well.indiana.attic.analyses.nahr;

import com.sun.net.httpserver.HttpServer;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.json.JSONObject;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.http.BasicHandler;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.*;
import java.net.InetSocketAddress;
import java.util.*;

public class MarkovianGraph extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (.fasta)")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="graph", shortName="g", doc="Graph (.ctx)")
    public HashSet<CortexGraph> GRAPHS;

    @Argument(fullName="port", shortName="p", doc="Port")
    public Integer PORT = 9001;

    @Output
    public PrintStream out;

    private Map<String, String> contigs = new HashMap<String, String>();

    private class PageHandler extends BasicHandler {
        @Override
        public String process(File page, Map<String, String> query) {
            try {
                if (page.exists()) {
                    InputStream is = new FileInputStream(page);

                    BufferedReader in = new BufferedReader(new InputStreamReader(is));
                    StringBuilder responseData = new StringBuilder();

                    String line;
                    while ((line = in.readLine()) != null) {
                        responseData.append(line).append("\n");
                    }

                    return responseData.toString();
                }
            } catch (IOException e) {
                return null;
            }

            return null;
        }
    }

    private class SearchHandler extends BasicHandler {
        @Override
        public String process(File page, Map<String, String> query) {
            if (query.containsKey("contigName") && contigs.containsKey(query.get("contigName"))) {
                String contig = contigs.get(query.get("contigName"));

                DirectedGraph<String, DefaultEdge> g = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);

                for (CortexGraph cg : GRAPHS) {
                    if (cg.getCortexFile().getName().contains("PG0063-C.ERR019060.infer.sorted.k41.ctx")) {
                        for (int i = 0; i <= contig.length() - cg.getKmerSize(); i++) {
                            String sk = contig.substring(i, i + cg.getKmerSize());

                            if (i == 0) {
                                for (int j = 0; j < cg.getKmerSize() - 1; j++) {
                                    String curId = sk.charAt(j) + "_" + j;
                                    String nextId = sk.charAt(j + 1) + "_" + (j+1);

                                    g.addVertex(curId);
                                    g.addVertex(nextId);
                                    g.addEdge(curId, nextId);
                                }
                            }

                            Set<String> prevKmers = CortexUtils.getPrevKmers(cg, sk);
                            for (String prevKmer : prevKmers) {
                                String prevId = prevKmer.charAt(0) + "_" + (i - 1);
                                String curId = sk.charAt(i) + "_" + i;

                                g.addVertex(prevId);
                                g.addVertex(curId);
                                g.addEdge(prevId, curId);
                            }

                            Set<String> nextKmers = CortexUtils.getNextKmers(cg, sk);
                            for (String nextKmer : nextKmers) {
                                String curId = sk.charAt(sk.length() - 1) + "_" + (i + cg.getKmerSize());
                                String nextId = nextKmer.charAt(nextKmer.length() - 1) + "_" + (i + cg.getKmerSize() + 1);

                                g.addVertex(curId);
                                g.addVertex(nextId);
                                g.addEdge(curId, nextId);
                            }
                        }
                    }
                }


                List<Map<String, String>> nodes = new ArrayList<Map<String, String>>();
                Map<String, Integer> nodePos = new HashMap<String, Integer>();
                int pos = 0;
                for (String node : g.vertexSet()) {
                    Map<String, String> entry = new HashMap<String, String>();
                    entry.put("name", node);
                    entry.put("group", "1");

                    nodes.add(entry);
                    nodePos.put(node, pos);
                    pos++;
                }

                List<Map<String, String>> links = new ArrayList<Map<String, String>>();
                for (DefaultEdge e : g.edgeSet()) {
                    String kmerSource = g.getEdgeSource(e);
                    String kmerTarget = g.getEdgeTarget(e);

                    int posSource = nodePos.get(kmerSource);
                    int posTarget = nodePos.get(kmerTarget);

                    Map<String, String> entry = new HashMap<String, String>();
                    //entry.put("source", String.valueOf(posSource));
                    //entry.put("target", String.valueOf(posTarget));
                    entry.put("source", kmerSource);
                    entry.put("target", kmerTarget);
                    entry.put("value", "1");
                    links.add(entry);
                }

                JSONObject jo = new JSONObject();
                jo.put("contig", contig);
                //jo.put("nodes", nodes);
                jo.put("links", links);

                return jo.toString();
            }

            return null;
        }
    }

    private void loadContigs() {
        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            contigs.put(rseq.getName(), new String(rseq.getBases()));
        }
    }

    @Override
    public void execute() {
        loadContigs();

        log.info("Loaded {} contigs", contigs.size());
        log.info("Loaded {} graphs", GRAPHS.size());
        for (CortexGraph cg : GRAPHS) {
            log.info("  {}: {}", cg.getCortexFile().getName(), cg.getKmerSize());
        }

        try {
            HttpServer server = HttpServer.create(new InetSocketAddress(PORT), 0);
            server.createContext("/", new PageHandler());
            server.createContext("/search", new SearchHandler());
            server.setExecutor(null);
            server.start();
        } catch (IOException e) {
            throw new IndianaException("Unable to start server", e);
        }
    }
}
