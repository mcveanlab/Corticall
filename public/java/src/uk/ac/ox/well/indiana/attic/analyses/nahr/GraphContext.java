package uk.ac.ox.well.indiana.attic.analyses.nahr;

import com.sun.net.httpserver.HttpServer;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.ivy.util.StringUtils;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.json.JSONObject;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.http.BasicHandler;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksMap;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.*;
import java.net.InetSocketAddress;
import java.util.*;

public class GraphContext extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (.fasta)")
    public File CONTIGS;

    @Argument(fullName="metrics", shortName="m", doc="Metrics (.metrics)")
    public File METRICS;

    @Argument(fullName="reference", shortName="r", doc="Reference (.fasta)")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="graph", shortName="g", doc="Graph (.ctx)")
    public LinkedHashSet<CortexGraph> GRAPHS;

    @Argument(fullName="links", shortName="l", doc="Links (.ctp)", required=false)
    public HashSet<CortexLinksMap> LINKS;

    @Argument(fullName="port", shortName="p", doc="Port")
    public Integer PORT = 9001;

    @Output
    public PrintStream out;

    private Map<String, String> contigs = new HashMap<String, String>();
    private Map<String, Map<String, String>> metrics = new HashMap<String, Map<String, String>>();

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
        private boolean missingFromGraph;
        private int cov;

        public CtxVertex(String kmer, int pos, VertexType vt, CortexRecord cr) {
            this.kmer = kmer;
            this.pos = pos + kmer.length() - 1;
            this.vertexType = vt;
            this.missingFromGraph = (cr == null);
            this.cov = (cr == null) ? 0 : cr.getCoverage(0);
        }

        public String getKmer() { return this.kmer; }
        public int getPos() { return this.pos; }
        public VertexType getVertexType() { return this.vertexType; }
        public String getBase() { return this.kmer.substring(kmer.length() - 1); }
        public boolean isMissingFromGraph() { return this.missingFromGraph; }
        public int getCoverage() { return this.cov; }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            CtxVertex ctxVertex = (CtxVertex) o;

            if (cov != ctxVertex.cov) return false;
            if (missingFromGraph != ctxVertex.missingFromGraph) return false;
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
            result = 31 * result + (missingFromGraph ? 1 : 0);
            result = 31 * result + cov;
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

    private Cigar cigarStringToCigar(String cs) {
        List<CigarElement> ces = new ArrayList<CigarElement>();

        int start = 0;
        int stop = 1;

        while (stop < cs.length()) {
            char c = cs.charAt(stop);

            if (c == 'M' || c == 'S' || c == 'H' || c == 'I' || c == 'D') {
                int length = Integer.valueOf(cs.substring(start, stop));
                CigarOperator co = CigarOperator.characterToEnum(c);

                ces.add(new CigarElement(length, co));

                start = stop + 1;
            }

            stop++;
        }

        return new Cigar(ces);
    }

    private class SearchHandler extends BasicHandler {
        @Override
        public String process(File page, Map<String, String> query) {
            loadContigs();

            if (query.get("contigName").matches("^[ACGT]+$")) {
                contigs.put("manual", query.get("contigName"));
                query.put("contigName", "manual");
            } else if (query.get("contigName").matches("^Pf3D7.+$")) {
                String[] pieces = query.get("contigName").split("[:-]");

                int start = Integer.valueOf(pieces[1].replaceAll(",", ""));
                int end = Integer.valueOf(pieces[2].replaceAll(",", ""));

                ReferenceSequence rseq = REF.getSubsequenceAt(pieces[0], start, end);
                contigs.put("manual", new String(rseq.getBases()));
                query.put("contigName", "manual");
            }

            if (query.containsKey("contigName") && contigs.containsKey(query.get("contigName")) && graphs.containsKey(query.get("graphName"))) {
                boolean showLinks = query.get("showLinks").equals("links_on");

                String contig = contigs.get(query.get("contigName"));
                String originalContig = contigs.get(query.get("contigName"));
                String refFormattedString = "";
                String kmerOrigin = "";

                if (metrics.containsKey(query.get("contigName"))) {
                    String[] loc = metrics.get(query.get("contigName")).get("canonicalLocus").split("[:-]");
                    if (!loc[0].equals("*")) {
                        boolean isRc = metrics.get(query.get("contigName")).get("isRcCanonical").equals("1");

                        if (isRc) {
                            contig = SequenceUtils.reverseComplement(contig);
                            originalContig = SequenceUtils.reverseComplement(originalContig);
                        }

                        int locStart = Integer.valueOf(loc[1]);
                        int locEnd = Integer.valueOf(loc[2]);

                        Cigar cigar = cigarStringToCigar(metrics.get(query.get("contigName")).get("cigarCanonical"));
                        if (cigar.getCigarElement(0).getOperator().equals(CigarOperator.S)) {
                            locStart -= cigar.getCigarElement(0).getLength();
                        }

                        if (cigar.getCigarElement(cigar.getCigarElements().size() - 1).getOperator().equals(CigarOperator.S)) {
                            locEnd += cigar.getCigarElement(cigar.getCigarElements().size() - 1).getLength();
                        }

                        String ref = new String(REF.getSubsequenceAt(loc[0], locStart, locEnd).getBases());

                        StringBuilder refFormatted = new StringBuilder();
                        int pos = 0;
                        for (CigarElement ce : cigar.getCigarElements()) {
                            CigarOperator co = ce.getOperator();
                            switch (co) {
                                case S:
                                    refFormatted.append(ref.substring(pos, pos + ce.getLength()).toLowerCase());
                                    break;
                                case M:
                                    refFormatted.append(ref.substring(pos, pos + ce.getLength()));
                                    break;
                                case I:
                                    refFormatted.append(StringUtils.repeat("-", ce.getLength()));
                                    break;
                            }

                            if (ce.getOperator().consumesReferenceBases()) {
                                pos += ce.getLength();
                            }
                        }

                        refFormattedString = refFormatted.toString();

                        kmerOrigin = metrics.get(query.get("contigName")).get("kmerOrigin");
                    }
                }

                CortexGraph cg = graphs.get(query.get("graphName"));

                String sampleName = cg.getColor(0).getSampleName();
                Set<CortexLinksMap> links = new HashSet<CortexLinksMap>();
                if (LINKS != null && !LINKS.isEmpty()) {
                    for (CortexLinksMap link : LINKS) {
                        if (sampleName.equals(link.getCortexLinks().getColor(0).getSampleName())) {
                            links.add(link);
                        }
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
                Set<String> usedPrevKmers = new HashSet<String>();
                usedPrevKmers.add(firstKmer);
                while (pks.size() == 1 && usedPrevKmers.size() <= 100) {
                    String kmer = pks.iterator().next();
                    firstFlank.insert(0, kmer.charAt(0));

                    if (usedPrevKmers.contains(kmer)) {
                        break;
                    }
                    usedPrevKmers.add(kmer);

                    pks = CortexUtils.getPrevKmers(cg, kmer);
                }

                StringBuilder lastFlank = new StringBuilder();
                String lastKmer = contig.substring(contig.length() - cg.getKmerSize(), contig.length());
                Set<String> nks = CortexUtils.getNextKmers(cg, lastKmer);
                Set<String> usedNextKmers = new HashSet<String>();
                usedNextKmers.add(lastKmer);
                while (nks.size() == 1 && usedNextKmers.size() <= 100) {
                    String kmer = nks.iterator().next();
                    lastFlank.append(kmer.charAt(kmer.length() - 1));

                    if (usedNextKmers.contains(kmer)) {
                        break;
                    }
                    usedNextKmers.add(kmer);

                    nks = CortexUtils.getNextKmers(cg, kmer);
                }

                contig = firstFlank.toString() + contig + lastFlank.toString();

                DirectedGraph<CtxVertex, MultiEdge> g = new DefaultDirectedGraph<CtxVertex, MultiEdge>(MultiEdge.class);
                for (int i = 0; i <= contig.length() - cg.getKmerSize(); i++) {
                    String curKmer = contig.substring(i, i + cg.getKmerSize());
                    CortexKmer ck = new CortexKmer(curKmer);
                    CtxVertex curVer = new CtxVertex(curKmer, i, contigKmers.contains(curKmer) ? VertexType.CONTIG : VertexType.CLIPPED, cg.findRecord(ck));

                    g.addVertex(curVer);

                    String expectedPrevKmer = (i > 0) ? contig.substring(i - 1, i - 1 + cg.getKmerSize()) : "";
                    String expectedNextKmer = (i < contig.length() - cg.getKmerSize()) ? contig.substring(i + 1, i + 1 + cg.getKmerSize()) : "";

                    Set<String> prevKmers = CortexUtils.getPrevKmers(cg, curKmer);
                    for (String prevKmer : prevKmers) {
                        if (!expectedPrevKmer.equals(prevKmer)) {
                            CortexKmer pk = new CortexKmer(prevKmer);
                            CtxVertex prevVer = new CtxVertex(prevKmer, i - 1, VertexType.IN, cg.findRecord(pk));

                            MultiEdge me = g.containsEdge(prevVer, curVer) ? g.getEdge(prevVer, curVer) : new MultiEdge();
                            me.addGraphName(cg.getCortexFile().getName());

                            g.addVertex(prevVer);
                            g.addEdge(prevVer, curVer, me);
                        }
                    }

                    Set<String> nextKmers = CortexUtils.getNextKmers(cg, curKmer);
                    for (String nextKmer : nextKmers) {
                        if (!expectedNextKmer.equals(nextKmer)) {
                            CortexKmer nk = new CortexKmer(nextKmer);
                            CtxVertex nextVer = new CtxVertex(nextKmer, i + 1, VertexType.OUT, cg.findRecord(nk));

                            MultiEdge me = g.containsEdge(curVer, nextVer) ? g.getEdge(curVer, nextVer) : new MultiEdge();
                            me.addGraphName(cg.getCortexFile().getName());

                            g.addVertex(nextVer);
                            g.addEdge(curVer, nextVer, me);
                        }
                    }
                }

                Set<Map<String, Object>> verticesWithLinks = new HashSet<Map<String, Object>>();
                DataFrame<String, String, Integer> hv = new DataFrame<String, String, Integer>(0);

                for (int q = 0; q <= contig.length() - cg.getKmerSize(); q++) {
                    //String sk = cv.getBinaryKmer();
                    String sk = contig.substring(q, q + cg.getKmerSize());
                    CortexKmer ck = new CortexKmer(sk);

                    for (CortexLinksMap link : links) {
                        if (link.containsKey(ck)) {
                            CortexLinksRecord clr = link.get(ck);
                            Map<String, Integer> lc = (!showLinks) ? new HashMap<String, Integer>() : CortexUtils.getKmersAndCoverageInLink(cg, sk, clr);

                            Map<String, Object> entry = new HashMap<String, Object>();
                            entry.put("kmer", sk);
                            entry.put("lc", lc);

                            verticesWithLinks.add(entry);

                            if (showLinks) {
                                for (CortexJunctionsRecord cjr : clr.getJunctions()) {
                                    List<String> lk = CortexUtils.getKmersInLink(cg, sk, cjr);

                                    for (int i = 0; i < lk.size(); i++) {
                                        String kili = lk.get(i);

                                        for (int j = 0; j < lk.size(); j++) {
                                            String kilj = lk.get(j);

                                            if (i != j) {
                                                hv.set(kili, kilj, hv.get(kili, kilj) + cjr.getCoverage(0));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                /*
                int hvMax = 0;
                Map<String, Integer> hvlin = new HashMap<String, Integer>();
                if (showLinks) {
                    for (String kili : hv.getRowNames()) {
                        for (String kilj : hv.getColNames()) {
                            int cov = hv.get(kili, kilj);

                            String id = kili + "_" + kilj;
                            hvlin.put(id, cov);

                            if (cov > hvMax) {
                                hvMax = cov;
                            }
                        }
                    }
                }
                */

                JSONObject jo = new JSONObject();
                jo.put("contig", contig);
                jo.put("originalContig", originalContig);
                jo.put("ref", refFormattedString);
                jo.put("kmerOrigin", kmerOrigin);
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
                    vm.put("missing", v.isMissingFromGraph());
                    vm.put("cov", v.getCoverage());

                    va.add(vm);
                }

                jo.put("vertices", va);
                jo.put("verticesWithLinks", verticesWithLinks);
                //jo.put("hvlin", hvlin);
                //jo.put("hvmax", hvMax);

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
                jo.put("cr", query.get("sk") + " not found");

                CortexGraph cg = graphs.get(query.get("graphName"));

                CortexKmer ck = new CortexKmer(query.get("sk"));

                CortexRecord cr = cg.findRecord(ck);
                if (cr != null) {
                    String text = cr.toString();
                    if (ck.isFlipped()) {
                        String info = query.get("sk");

                        for (int coverage : cr.getCoverages()) {
                            info += " " + coverage;
                        }

                        for (String edge : cr.getEdgeAsStrings()) {
                            info += " " + SequenceUtils.reverseComplement(edge);
                        }

                        text = info;
                    }

                    String sampleName = cg.getColor(0).getSampleName();
                    if (LINKS != null && !LINKS.isEmpty()) {
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

                                clr = CortexUtils.orientLinksRecord(query.get("sk"), clr);

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
        if (contigs == null || contigs.size() == 0) {
            contigs = new HashMap<String, String>();

            TableReader tr = new TableReader(METRICS);
            metrics = new HashMap<String, Map<String, String>>();

            for (Map<String, String> te : tr) {
                String contigName = te.get("contigName");
                String seq = te.get("seq");

                contigs.put(contigName, seq);

                metrics.put(contigName, te);
            }
        }

        FastaSequenceFile fasta = new FastaSequenceFile(CONTIGS, false);
        ReferenceSequence rseq;
        while ((rseq = fasta.nextSequence()) != null) {
            String[] name = rseq.getName().split("\\s+");
            contigs.put(name[0], new String(rseq.getBases()));
        }
    }

    private Map<String, CortexGraph> graphs = new HashMap<String, CortexGraph>();

    @Override
    public void execute() {
        Map<String, Set<CortexLinksMap>> links = new HashMap<String, Set<CortexLinksMap>>();
        Map<String, Set<String>> linkFilenames = new HashMap<String, Set<String>>();
        if (LINKS != null && !LINKS.isEmpty()) {
            for (CortexLinksMap cl : LINKS) {
                String name = cl.getCortexLinks().getColor(0).getSampleName();

                if (!links.containsKey(name)) {
                    links.put(name, new HashSet<CortexLinksMap>());
                    linkFilenames.put(name, new HashSet<String>());
                }

                links.get(name).add(cl);
                linkFilenames.get(name).add(cl.getCortexLinks().getCortexLinksFile().getAbsolutePath());
            }
            log.info("Loaded {} links", LINKS.size());

            for (CortexLinksMap clm : LINKS) {
                log.info("  {}: {}", clm.getCortexLinks().getCortexLinksFile().getAbsolutePath(), clm.getCortexLinks().getNumLinks());
            }
        }

        log.info("Loaded {} graphs", GRAPHS.size());
        for (CortexGraph cg : GRAPHS) {
            String name = cg.getColor(0).getSampleName();
            graphs.put(cg.getCortexFile().getName(), cg);

            log.info("  {}: {} {}", cg.getCortexFile().getName(), cg.getKmerSize(), linkFilenames.get(name));
        }

        loadContigs();
        log.info("Loaded {} contigs", contigs.size());

        log.info("Listening on port {}", PORT);

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
