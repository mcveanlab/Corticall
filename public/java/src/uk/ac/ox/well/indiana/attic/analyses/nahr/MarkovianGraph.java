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
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksMap;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.*;
import java.net.InetSocketAddress;
import java.util.*;

public class MarkovianGraph extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (.fasta)")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="graph", shortName="g", doc="Graph (.ctx)")
    public LinkedHashSet<CortexGraph> GRAPHS;

    @Argument(fullName="links", shortName="l", doc="Links (.ctp)")
    public HashSet<CortexLinksMap> LINKS;

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

    private class GraphList extends BasicHandler {
        @Override
        public String process(File page, Map<String, String> query) {
            Set<String> graphList = new TreeSet<String>();

            for (CortexGraph cg : GRAPHS) {
                graphList.add(cg.getCortexFile().getName());
            }

            JSONObject jo = new JSONObject();
            jo.put("graphList", graphList);

            return jo.toString();
        }
    }

    private class MultiEdge extends DefaultEdge {
        private Set<String> graphNames = new HashSet<String>();

        public void addGraphName(String graphName) {
            graphNames.add(graphName);
        }

        public Set<String> getGraphNames() {
            return graphNames;
        }
    }

    private class SearchHandler extends BasicHandler {
        @Override
        public String process(File page, Map<String, String> query) {
            if (query.containsKey("contigName") && contigs.containsKey(query.get("contigName")) && graphs.containsKey(query.get("graphName"))) {
                String contig = contigs.get(query.get("contigName"));
                CortexGraph cg = graphs.get(query.get("graphName"));

                String sampleName = cg.getColor(0).getSampleName();
                Set<CortexLinksMap> links = new HashSet<CortexLinksMap>();
                for (CortexLinksMap link : LINKS) {
                    if (sampleName.equals(link.getCortexLinks().getColor(0).getSampleName())) {
                        links.add(link);
                    }
                }

                DirectedGraph<String, MultiEdge> g = new DefaultDirectedGraph<String, MultiEdge>(MultiEdge.class);
                Set<String> nodesWithLinks = new HashSet<String>();

                Set<List<String>> kmersInLinks = new HashSet<List<String>>();

                for (int i = 0; i <= contig.length() - cg.getKmerSize(); i++) {
                    String curKmer = contig.substring(i, i + cg.getKmerSize());
                    CortexKmer ck = new CortexKmer(curKmer);

                    Set<String> prevKmers = CortexUtils.getPrevKmers(cg, curKmer);
                    for (String prevKmer : prevKmers) {
                        g.addVertex(prevKmer);
                        g.addVertex(curKmer);

                        MultiEdge me = g.containsEdge(prevKmer, curKmer) ? g.getEdge(prevKmer, curKmer) : new MultiEdge();
                        me.addGraphName(cg.getCortexFile().getName());

                        g.addEdge(prevKmer, curKmer, me);
                    }

                    Set<String> nextKmers = CortexUtils.getNextKmers(cg, curKmer);
                    for (String nextKmer : nextKmers) {
                        g.addVertex(curKmer);
                        g.addVertex(nextKmer);

                        MultiEdge me = g.containsEdge(curKmer, nextKmer) ? g.getEdge(curKmer, nextKmer) : new MultiEdge();
                        me.addGraphName(cg.getCortexFile().getName());

                        g.addEdge(curKmer, nextKmer, me);
                    }

                    for (CortexLinksMap link : links) {
                        if (link.containsKey(ck)) {
                            nodesWithLinks.add(curKmer);

                            log.info("link: {} {} {}", cg.getKmerSize(), i, curKmer);

                            for (CortexJunctionsRecord cjr : link.get(ck).getJunctions()) {
                                kmersInLinks.add(CortexUtils.getKmersInLink(cg, curKmer, cjr));
                            }
                        }
                    }
                }

                List<Map<String, Object>> gl = new ArrayList<Map<String, Object>>();
                for (MultiEdge e : g.edgeSet()) {
                    String kmerSource = g.getEdgeSource(e);
                    String kmerTarget = g.getEdgeTarget(e);

                    Map<String, Object> entry = new HashMap<String, Object>();
                    entry.put("source", kmerSource);
                    entry.put("target", kmerTarget);
                    entry.put("graphs", e.getGraphNames());
                    entry.put("value", e.getGraphNames().size());
                    gl.add(entry);
                }

                JSONObject jo = new JSONObject();
                jo.put("contig", contig);
                jo.put("links", gl);
                jo.put("nodesWithLinks", nodesWithLinks);
                jo.put("kmersInLinks", kmersInLinks);

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

    private Map<String, CortexGraph> graphs = new HashMap<String, CortexGraph>();

    @Override
    public void execute() {
        Map<String, Set<CortexLinksMap>> links = new HashMap<String, Set<CortexLinksMap>>();
        Map<String, Set<String>> linkFilenames = new HashMap<String, Set<String>>();
        for (CortexLinksMap cl : LINKS) {
            String name = cl.getCortexLinks().getColor(0).getSampleName();

            if (!links.containsKey(name)) {
                links.put(name, new HashSet<CortexLinksMap>());
                linkFilenames.put(name, new HashSet<String>());
            }

            links.get(name).add(cl);
            linkFilenames.get(name).add(cl.getCortexLinks().getCortexLinksFile().getName());
        }
        log.info("Loaded {} links", LINKS.size());

        log.info("Loaded {} graphs", GRAPHS.size());
        for (CortexGraph cg : GRAPHS) {
            String name = cg.getColor(0).getSampleName();
            graphs.put(cg.getCortexFile().getName(), cg);

            log.info("  {}: {} {}", cg.getCortexFile().getName(), cg.getKmerSize(), linkFilenames.get(name));
        }

        loadContigs();
        log.info("Loaded {} contigs", contigs.size());

        try {
            HttpServer server = HttpServer.create(new InetSocketAddress(PORT), 0);
            server.createContext("/", new PageHandler());
            server.createContext("/graphlist", new GraphList());
            server.createContext("/search", new SearchHandler());
            server.setExecutor(null);
            server.start();
        } catch (IOException e) {
            throw new IndianaException("Unable to start server", e);
        }
    }
}
