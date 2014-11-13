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
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.*;
import java.net.InetSocketAddress;
import java.util.*;

public class MarkovianGraph extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (.fasta)")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="graph", shortName="g", doc="Graph (.ctx)")
    public LinkedHashSet<CortexGraph> GRAPHS;

    @Argument(fullName="links", shortName="l", doc="Links (.ctp)")
    public HashSet<CortexLinks> LINKS;

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
            if (query.containsKey("contigName") && contigs.containsKey(query.get("contigName"))) {
                String contig = contigs.get(query.get("contigName"));

                DirectedGraph<String, MultiEdge> g = new DefaultDirectedGraph<String, MultiEdge>(MultiEdge.class);
                Set<String> nodesWithLinks = new HashSet<String>();

                for (CortexGraph cg : GRAPHS) {
                    String sampleName = cg.getColor(0).getSampleName();
                    Set<CortexLinks> links = new HashSet<CortexLinks>();
                    Set<String> linkNames = new HashSet<String>();
                    for (CortexLinks link : LINKS) {
                        if (sampleName.equals(link.getColor(0).getSampleName())) {
                            links.add(link);
                            linkNames.add(link.getCortexLinksFile().getName());
                        }
                    }

                    for (int i = 0; i <= contig.length() - cg.getKmerSize(); i++) {
                        String sk = contig.substring(i, i + cg.getKmerSize());
                        CortexKmer ck = new CortexKmer(sk);

                        for (CortexLinks link : links) {
                            if (linkRecords.get(link).containsKey(ck)) {
                                String linkId = sk.charAt(sk.length() - 1) + "_" + (i + cg.getKmerSize() - 1);
                                nodesWithLinks.add(linkId);

                                log.info("link: {} {} {}", cg.getKmerSize(), i, linkId);
                            }
                        }

                        if (i == 0) {
                            for (int j = 0; j < cg.getKmerSize() - 1; j++) {
                                String curId = sk.charAt(j) + "_" + j;
                                String nextId = sk.charAt(j + 1) + "_" + (j + 1);

                                g.addVertex(curId);
                                g.addVertex(nextId);

                                MultiEdge me = g.containsEdge(curId, nextId) ? g.getEdge(curId, nextId) : new MultiEdge();
                                me.addGraphName(cg.getCortexFile().getName());

                                g.addEdge(curId, nextId, me);
                            }
                        }

                        Set<String> prevKmers = CortexUtils.getPrevKmers(cg, sk);
                        for (String prevKmer : prevKmers) {
                            String prevId = prevKmer.charAt(0) + "_" + (i - 1);
                            String curId = sk.charAt(0) + "_" + i;

                            g.addVertex(prevId);
                            g.addVertex(curId);

                            MultiEdge me = g.containsEdge(prevId, curId) ? g.getEdge(prevId, curId) : new MultiEdge();
                            me.addGraphName(cg.getCortexFile().getName());

                            g.addEdge(prevId, curId, me);
                        }

                        Set<String> nextKmers = CortexUtils.getNextKmers(cg, sk);
                        for (String nextKmer : nextKmers) {
                            String curId = sk.charAt(sk.length() - 1) + "_" + (i + cg.getKmerSize() - 1);
                            String nextId = nextKmer.charAt(nextKmer.length() - 1) + "_" + (i + cg.getKmerSize());

                            g.addVertex(curId);
                            g.addVertex(nextId);

                            MultiEdge me = g.containsEdge(curId, nextId) ? g.getEdge(curId, nextId) : new MultiEdge();
                            me.addGraphName(cg.getCortexFile().getName());

                            g.addEdge(curId, nextId, me);
                        }
                    }
                }

                Set<String> contigNodes = new HashSet<String>();
                for (int i = 0; i < contig.length(); i++) {
                    String curId = contig.charAt(i) + "_" + String.valueOf(i);

                    contigNodes.add(curId);
                }

                List<Map<String, Object>> links = new ArrayList<Map<String, Object>>();
                for (MultiEdge e : g.edgeSet()) {
                    String kmerSource = g.getEdgeSource(e);
                    String kmerTarget = g.getEdgeTarget(e);

                    Map<String, Object> entry = new HashMap<String, Object>();
                    entry.put("source", kmerSource);
                    entry.put("target", kmerTarget);
                    entry.put("graphs", e.getGraphNames());
                    entry.put("fixed", contigNodes.contains(kmerSource) && contigNodes.contains(kmerTarget) ? "true" : "false");
                    entry.put("value", e.getGraphNames().size());
                    links.add(entry);
                }

                JSONObject jo = new JSONObject();
                jo.put("contig", contig);
                jo.put("links", links);
                jo.put("nodesWithLinks", nodesWithLinks);

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

    private Map<CortexLinks, Map<CortexKmer, CortexLinksRecord>> linkRecords = new HashMap<CortexLinks, Map<CortexKmer, CortexLinksRecord>>();

    private void loadLinkRecords(CortexLinks cl) {
        if (!linkRecords.containsKey(cl)) {
            linkRecords.put(cl, new HashMap<CortexKmer, CortexLinksRecord>());

            for (CortexLinksRecord clr : cl) {
                linkRecords.get(cl).put(clr.getKmer(), clr);
            }
        }
    }

    @Override
    public void execute() {
        loadContigs();
        log.info("Loaded {} contigs", contigs.size());

        Map<String, Set<CortexLinks>> links = new HashMap<String, Set<CortexLinks>>();
        for (CortexLinks cl : LINKS) {
            String name = cl.getColor(0).getSampleName();

            if (!links.containsKey(name)) {
                links.put(name, new HashSet<CortexLinks>());
            }

            links.get(name).add(cl);

            loadLinkRecords(cl);
        }

        log.info("Loaded {} links", linkRecords.size());

        log.info("Loaded {} graphs", GRAPHS.size());
        for (CortexGraph cg : GRAPHS) {
            String name = cg.getColor(0).getSampleName();

            Set<String> linkFiles = new HashSet<String>();
            if (links.containsKey(name)) {
                for (CortexLinks cl : links.get(name)) {
                    linkFiles.add(cl.getCortexLinksFile().getName());
                }
            }

            log.info("  {}: {} {}", cg.getCortexFile().getName(), cg.getKmerSize(), linkFiles);
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
