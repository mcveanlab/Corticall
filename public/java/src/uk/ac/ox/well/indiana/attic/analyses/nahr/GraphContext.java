package uk.ac.ox.well.indiana.attic.analyses.nahr;

import com.sun.net.httpserver.HttpServer;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.json.JSONArray;
import org.json.JSONObject;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.http.BasicHandler;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksMap;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.*;
import java.net.InetSocketAddress;
import java.util.*;

public class GraphContext extends Module {
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

    public enum VertexType { CONTIG, OUT, IN, CLIPPED }

    private class CtxVertex {
        private String kmer;
        private int pos;
        private VertexType vertexType;

        public CtxVertex(String kmer, int pos, VertexType vt) {
            this.kmer = kmer;
            this.pos = pos + kmer.length() - 1;
            this.vertexType = vt;
        }

        public String getKmer() { return this.kmer; }
        public int getPos() { return this.pos; }
        public VertexType getVertexType() { return this.vertexType; }
        public String getBase() { return this.kmer.substring(kmer.length() - 1); }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            CtxVertex ctxVertex = (CtxVertex) o;

            if (pos != ctxVertex.pos) return false;
            if (kmer != null ? !kmer.equals(ctxVertex.kmer) : ctxVertex.kmer != null) return false;
            if (vertexType != ctxVertex.vertexType) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = kmer != null ? kmer.hashCode() : 0;
            result = 31 * result + pos;
            result = 31 * result + (vertexType != null ? vertexType.hashCode() : 0);
            return result;
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

                Set<String> contigKmers = new HashSet<String>();
                for (int i = 0; i <= contig.length() - cg.getKmerSize(); i++) {
                    String curKmer = contig.substring(i, i + cg.getKmerSize());

                    contigKmers.add(curKmer);
                }

                StringBuilder firstFlank = new StringBuilder();
                String firstKmer = contig.substring(0, cg.getKmerSize());
                Set<String> pks = CortexUtils.getPrevKmers(cg, firstKmer);
                while (pks.size() == 1) {
                    String kmer = pks.iterator().next();
                    firstFlank.insert(0, kmer.charAt(0));

                    pks = CortexUtils.getPrevKmers(cg, kmer);
                }

                StringBuilder lastFlank = new StringBuilder();
                String lastKmer = contig.substring(contig.length() - cg.getKmerSize(), contig.length());
                Set<String> nks = CortexUtils.getNextKmers(cg, lastKmer);
                while (nks.size() == 1) {
                    String kmer = nks.iterator().next();
                    lastFlank.append(kmer.charAt(kmer.length() - 1));

                    nks = CortexUtils.getNextKmers(cg, kmer);
                }

                contig = firstFlank.toString() + contig + lastFlank.toString();

                DirectedGraph<CtxVertex, MultiEdge> g = new DefaultDirectedGraph<CtxVertex, MultiEdge>(MultiEdge.class);
                for (int i = 0; i <= contig.length() - cg.getKmerSize(); i++) {
                    String curKmer = contig.substring(i, i + cg.getKmerSize());
                    CtxVertex curVer = new CtxVertex(curKmer, i, contigKmers.contains(curKmer) ? VertexType.CONTIG : VertexType.CLIPPED);

                    g.addVertex(curVer);

                    String expectedPrevKmer = (i > 0) ? contig.substring(i - 1, i - 1 + cg.getKmerSize()) : "";
                    String expectedNextKmer = (i < contig.length() - cg.getKmerSize()) ? contig.substring(i + 1, i + 1 + cg.getKmerSize()) : "";

                    Set<String> prevKmers = CortexUtils.getPrevKmers(cg, curKmer);
                    for (String prevKmer : prevKmers) {
                        if (!expectedPrevKmer.equals(prevKmer)) {
                            CtxVertex prevVer = new CtxVertex(prevKmer, i - 1, VertexType.IN);

                            MultiEdge me = g.containsEdge(prevVer, curVer) ? g.getEdge(prevVer, curVer) : new MultiEdge();
                            me.addGraphName(cg.getCortexFile().getName());

                            g.addVertex(prevVer);
                            g.addEdge(prevVer, curVer, me);
                        }
                    }

                    Set<String> nextKmers = CortexUtils.getNextKmers(cg, curKmer);
                    for (String nextKmer : nextKmers) {
                        if (!expectedNextKmer.equals(nextKmer)) {
                            CtxVertex nextVer = new CtxVertex(nextKmer, i + 1, VertexType.OUT);

                            MultiEdge me = g.containsEdge(curVer, nextVer) ? g.getEdge(curVer, nextVer) : new MultiEdge();
                            me.addGraphName(cg.getCortexFile().getName());

                            g.addVertex(nextVer);
                            g.addEdge(curVer, nextVer, me);
                        }
                    }
                }

                Set<Map<String, Object>> verticesWithLinks = new HashSet<Map<String, Object>>();
                Set<Map<String, Object>> kmersInLinks = new HashSet<Map<String, Object>>();

                for (CtxVertex cv : g.vertexSet()) {
                    String sk = cv.getKmer();
                    CortexKmer ck = new CortexKmer(sk);
                    for (CortexLinksMap link : links) {
                        if (link.containsKey(ck)) {
                            Map<String, Object> entry = new HashMap<String, Object>();
                            entry.put("kmer", sk);
                            entry.put("flipped", ck.isFlipped());
                            verticesWithLinks.add(entry);

                            CortexLinksRecord clr = link.get(ck);
                            for (CortexJunctionsRecord cjr : clr.getJunctions()) {
                                int cov = cjr.getCoverage(0);
                                List<String> kmersInLink = CortexUtils.getKmersInLink(cg, sk, cjr);

                                Map<String, Object> kl = new HashMap<String, Object>();
                                kl.put("kmers", kmersInLink);
                                kl.put("cov", cov);

                                kmersInLinks.add(kl);
                            }
                        }
                    }
                }

                JSONObject jo = new JSONObject();
                jo.put("contig", contig);
                jo.put("kmerSize", cg.getKmerSize());
                jo.put("clipStart", firstFlank.length());
                jo.put("clipEnd", contig.length() - lastFlank.length());

                List<Map<String, Object>> va = new ArrayList<Map<String, Object>>();
                for (CtxVertex v : g.vertexSet()) {
                    Map<String, Object> vm = new HashMap<String, Object>();
                    vm.put("base", v.getBase());
                    vm.put("kmer", v.getKmer());
                    vm.put("pos", v.getPos());
                    vm.put("type", v.getVertexType().name());

                    va.add(vm);
                }

                jo.put("vertices", va);
                jo.put("verticesWithLinks", verticesWithLinks);
                jo.put("kmersInLinks", kmersInLinks);

                return jo.toString();
            }

            return null;
        }
    }

    private class CtxRecordHandler extends BasicHandler {
        @Override
        public String process(File page, Map<String, String> query) {
            JSONObject jo = new JSONObject();
            jo.put("cr", "not found");

            if (graphs.containsKey(query.get("graphName")) && query.containsKey("sk")) {
                CortexGraph cg = graphs.get(query.get("graphName"));

                CortexKmer ck = new CortexKmer(query.get("sk"));

                CortexRecord cr = cg.findRecord(ck);
                if (cr != null) {
                    String text = cr.toString();

                    String sampleName = cg.getColor(0).getSampleName();
                    for (CortexLinksMap link : LINKS) {
                        if (sampleName.equals(link.getCortexLinks().getColor(0).getSampleName())) {
                            if (link.containsKey(ck)) {
                                CortexLinksRecord clr = link.get(ck);
                                int cov = 0;
                                for (CortexJunctionsRecord cjr : clr.getJunctions()) {
                                    cov += cjr.getCoverage(0);
                                }

                                text += " (" + cov + " links)";
                            }
                        }
                    }

                    jo.put("cr", text);
                }

            }

            return jo.toString();
        }
    }

    private class LinksHandler extends BasicHandler {
        @Override
        public String process(File page, Map<String, String> query) {
            JSONObject jo = new JSONObject();
            jo.put("tip", "");

            if (graphs.containsKey(query.get("graphName")) && query.containsKey("sk")) {
                CortexGraph cg = graphs.get(query.get("graphName"));

                CortexKmer ck = new CortexKmer(query.get("sk"));

                CortexRecord cr = cg.findRecord(ck);
                if (cr != null) {
                    String sampleName = cg.getColor(0).getSampleName();
                    for (CortexLinksMap link : LINKS) {
                        if (sampleName.equals(link.getCortexLinks().getColor(0).getSampleName())) {
                            if (link.containsKey(ck)) {
                                CortexLinksRecord clr = link.get(ck);

                                if (ck.isFlipped()) {
                                    clr = CortexUtils.flipLinksRecord(clr);
                                }

                                jo.put("tip", clr.toString().replaceAll("\n", ", "));

                                break;
                            }
                        }
                    }
                }

            }

            return jo.toString();
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
            server.createContext("/cr", new CtxRecordHandler());
            server.createContext("/links", new LinksHandler());
            server.setExecutor(null);
            server.start();
        } catch (IOException e) {
            throw new IndianaException("Unable to start server", e);
        }
    }
}
