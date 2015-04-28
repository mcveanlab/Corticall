package uk.ac.ox.well.indiana.commands.cortex;

import com.sun.net.httpserver.Headers;
import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlContext;
import org.apache.commons.jexl2.JexlEngine;
import org.apache.commons.jexl2.MapContext;
import org.json.JSONArray;
import org.json.JSONObject;
import org.yaml.snakeyaml.Yaml;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.*;
import java.net.InetSocketAddress;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

@Description(text="Starts the server for visualizing assembly data")
public class server extends Module {
    @Argument(fullName="manifest", shortName="m", doc="Manifest file (YAML)")
    public File MANIFEST;

    @Argument(fullName="port", shortName="p", doc="Port")
    public Integer PORT = 9000;

    private class Manifest {
        private Map<String, Object> manifestMap;

        public Manifest(File manifestFile) {
            Yaml yaml = new Yaml();
            try {
                manifestMap = (Map<String, Object>) yaml.load(new FileInputStream(manifestFile));
            } catch (FileNotFoundException e) {
                throw new IndianaException("Could not parse server configuration file, '" + manifestFile.getAbsolutePath() + "'", e);
            }
        }

        public Set<String> getCrosses() { return manifestMap.keySet(); }

        public Set<String> getSamples(String cross) {
            Map sampleMap = (Map) manifestMap.get(cross);
            return new TreeSet<String>(sampleMap.keySet());
        }

        public Set<String> getAccessions(String cross, String sample) {
            Map sampleMap = (Map) manifestMap.get(cross);
            Map accMap = (Map) sampleMap.get(sample);
            return new TreeSet<String>(accMap.keySet());
        }

        private Map getItems(String cross, String sample, String accession) {
            Map sampleMap = (Map) manifestMap.get(cross);
            Map accMap = (Map) sampleMap.get(sample);
            return (Map) accMap.get(accession);
        }

        public Map<String, SAMFileReader> getBams(String cross, String sample, String accession) {
            Map bamMap = (Map) getItems(cross, sample, accession).get("bams");

            Map<String, SAMFileReader> bamFileMap = new HashMap<String, SAMFileReader>();
            for (Object key : bamMap.keySet()) {
                File file = new File((String) bamMap.get(key));
                SAMFileReader sfr = new SAMFileReader(file);
                bamFileMap.put((String) key, sfr);
            }

            return bamFileMap;
        }

        public Map<String, GFF3> getGffs(String cross, String sample, String accession) {
            Map gffMap = (Map) getItems(cross, sample, accession).get("gffs");

            Map<String, GFF3> gffFileMap = new HashMap<String, GFF3>();
            for (Object key : gffMap.keySet()) {
                File file = new File((String) gffMap.get(key));
                gffFileMap.put((String) key, new GFF3(file));
            }

            return gffFileMap;
        }

        public Map<String, IndexedFastaSequenceFile> getRefs(String cross, String sample, String accession) {
            Map refMap = (Map) getItems(cross, sample, accession).get("refs");
            Map<String, IndexedFastaSequenceFile> refFileMap = new HashMap<String, IndexedFastaSequenceFile>();

            try {
                for (Object key : refMap.keySet()) {
                    File file = new File((String) refMap.get(key));
                    refFileMap.put((String) key, new IndexedFastaSequenceFile(file));
                }
            } catch (IOException e) {
                throw new IndianaException("Unable to load fasta", e);
            }

            return refFileMap;
        }

        public Set<String> getGraphKmerSizes(String cross, String sample, String accession) {
            Map graphsMap = (Map) getItems(cross, sample, accession).get("graphs");

            return new TreeSet<String>(graphsMap.keySet());
        }

        public Set<String> getGraphLabels(String cross, String sample, String accession, String kmerSize) {
            Map graphsMap = (Map) getItems(cross, sample, accession).get("graphs");
            Map graphMap = (Map) graphsMap.get(kmerSize);

            Set<String> graphLabelMap = new TreeSet<String>(graphMap.keySet());

            return graphLabelMap;
        }

        public Map<String, File> getGraphs(String cross, String sample, String accession, String kmerSize, String label) {
            Map graphsMap = (Map) ((Map) ((Map) getItems(cross, sample, accession).get("graphs")).get(kmerSize)).get(label);

            Map<String, File> graphs = new HashMap<String, File>();
            for (Object key : graphsMap.keySet()) {
                String skey = (String) key;

                if (skey.equals("graph")) {
                    File graph = new File((String) graphsMap.get(skey));

                    graphs.put(skey, graph);
                }
            }

            return graphs;
        }

        public Map<String, File> getLinks(String cross, String sample, String accession, String kmerSize, String label) {
            Map graphsMap = (Map) ((Map) ((Map) getItems(cross, sample, accession).get("graphs")).get(kmerSize)).get(label);

            Map<String, File> links = new HashMap<String, File>();
            for (Object key : graphsMap.keySet()) {
                String skey = (String) key;

                if (skey.startsWith("links_")) {
                    File link = new File((String) graphsMap.get(skey));

                    links.put(skey, link);
                }
            }

            return links;
        }

        public File getContigs(String cross, String sample, String accession) {
            Map contigsMap = (Map) getItems(cross, sample, accession).get("contigs");

            File contigs = null;

            for (Object key : contigsMap.keySet()) {
                String skey = (String) key;
                contigs = new File((String) contigsMap.get(skey));
            }

            return contigs;
        }

        public File getMetrics(String cross, String sample, String accession) {
            Map metricsMap = (Map) getItems(cross, sample, accession).get("metrics");

            File metrics = null;

            for (Object key : metricsMap.keySet()) {
                String skey = (String) key;
                //metrics = new File((String) metricsMap.get(skey));

                Map metricMap = (Map) metricsMap.get(skey);
                metrics = new File((String) metricMap.get("metrics"));
            }

            return metrics;
        }

        public File getCircos(String cross, String sample, String accession) {
            Map circosMap = (Map) getItems(cross, sample, accession).get("circos");

            File circos = null;

            for (Object key : circosMap.keySet()) {
                String skey = (String) key;
                //circos = new File((String) circosMap.get(skey));

                Map metricMap = (Map) circosMap.get(skey);
                circos = new File((String) metricMap.get("circos"));
            }

            return circos;
        }

        public int numSamples(String cross) {
            return getSamples(cross).size();
        }

        public int numAccessions(String cross) {
            Set<String> accessions = new HashSet<String>();

            for (String sample : getSamples(cross)) {
                accessions.addAll(getAccessions(cross, sample));
            }

            return accessions.size();
        }

        public int numKmerSizes(String cross) {
            Set<String> kmerSizes = new HashSet<String>();

            for (String sample : getSamples(cross)) {
                for (String accession : getAccessions(cross, sample)) {
                    for (String kmerSize : getGraphKmerSizes(cross, sample, accession)) {
                        kmerSizes.add(kmerSize);
                    }
                }
            }

            return kmerSizes.size();
        }

        public int numGraphs(String cross) {
            Set<File> graphs = new HashSet<File>();

            for (String sample : getSamples(cross)) {
                for (String accession : getAccessions(cross, sample)) {
                    for (String kmerSize : getGraphKmerSizes(cross, sample, accession)) {
                        for (String label : getGraphLabels(cross, sample, accession, kmerSize)) {
                            Map<String, File> graphFiles = getGraphs(cross, sample, accession, kmerSize, label);

                            graphs.addAll(graphFiles.values());
                        }
                    }
                }
            }

            return graphs.size();
        }

        public int numLinks(String cross) {
            Set<File> links = new HashSet<File>();

            for (String sample : getSamples(cross)) {
                for (String accession : getAccessions(cross, sample)) {
                    for (String kmerSize : getGraphKmerSizes(cross, sample, accession)) {
                        for (String label : getGraphLabels(cross, sample, accession, kmerSize)) {
                            Map<String, File> linkFiles = getLinks(cross, sample, accession, kmerSize, label);

                            links.addAll(linkFiles.values());
                        }
                    }
                }
            }

            return links.size();
        }
    }

    private Manifest m;

    private JexlEngine initializeJexlEngine() {
        JexlEngine je = new JexlEngine();
        je.setCache(512);
        je.setLenient(false);
        je.setSilent(false);

        return je;
    }

    private JexlEngine je;

    private abstract class BaseHandler implements HttpHandler {
        public Map<String, String> query(String query) {
            Map<String, String> result = new HashMap<String, String>();

            query = query.replaceAll("&&", "<and>");

            for (String param : query.split("&")) {
                String pair[] = param.split("=", 2);
                if (pair.length > 1) {
                    result.put(pair[0], pair[1].replaceAll("<and>", "&&"));
                } else {
                    result.put(pair[0], "");
                }
            }

            return result;
        }

        public void write(HttpExchange httpExchange, String response) throws IOException {
            write(httpExchange, 200, response);
        }

        public void write(HttpExchange httpExchange, int code, String response) throws IOException {
            httpExchange.sendResponseHeaders(code, response.length());

            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();
        }
    }

    private class PageHandler extends BaseHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            String page = "/html/" + httpExchange.getRequestURI();
            File f = new File("./public" + page);

            int code = 200;
            String response;

            if (f.exists()) {
                InputStream is;

                if (log.isDebugEnabled()) {
                    is = new FileInputStream(f);
                } else {
                    is = this.getClass().getResourceAsStream(page);
                }

                BufferedReader in = new BufferedReader(new InputStreamReader(is));
                StringBuilder responseData = new StringBuilder();

                String line;
                while ((line = in.readLine()) != null) {
                    responseData.append(line).append("\n");
                }

                Headers h = httpExchange.getResponseHeaders();
                if      (page.endsWith(".js"))   { h.add("Content-Type", "text/javascript"); }
                else if (page.endsWith(".css"))  { h.add("Content-Type", "text/css");        }
                else if (page.endsWith(".html")) { h.add("Content-Type", "text/html");       }

                response = responseData.toString();
            } else {
                code = 404;
                response = "Page not found.";
            }

            write(httpExchange, code, response);
        }
    }

    private class CrossesList extends BaseHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Set<String> crosses = m.getCrosses();

            JSONObject jo = new JSONObject();
            jo.put("crosses", crosses);

            write(httpExchange, jo.toString());
        }
    }

    private class SamplesList extends BaseHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = query(httpExchange.getRequestURI().getQuery());

            String cross = query.get("cross");

            Set<String> samples = m.getSamples(cross);

            JSONObject jo = new JSONObject();
            jo.put("samples", samples);

            write(httpExchange, jo.toString());
        }
    }

    private class AccessionsList extends BaseHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = query(httpExchange.getRequestURI().getQuery());

            String cross = query.get("cross");
            String sample = query.get("sample");

            Set<String> accessions = m.getAccessions(cross, sample);

            JSONObject jo = new JSONObject();
            jo.put("accessions", accessions);

            write(httpExchange, jo.toString());
        }
    }

    private class NumContigs extends BaseHandler {
        private Map<File, Integer> numContigsCache = new HashMap<File, Integer>();

        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = query(httpExchange.getRequestURI().getQuery());

            String cross = query.get("cross");
            String sample = query.get("sample");
            String accession = query.get("accession");

            File metricsFile = m.getMetrics(cross, sample, accession);

            if (!numContigsCache.containsKey(metricsFile)) {
                TableReader tr = new TableReader(metricsFile);

                int numContigs = 0;
                for (Map<String, String> te : tr) {
                    numContigs++;
                }

                numContigsCache.put(metricsFile, numContigs);
            }

            JSONObject jo = new JSONObject();
            jo.put("numContigs", numContigsCache.get(metricsFile));

            write(httpExchange, jo.toString());
        }
    }

    private class Circos extends BaseHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = query(httpExchange.getRequestURI().getQuery());

            String cross = query.get("cross");
            String sample = query.get("sample");
            String accession = query.get("accession");

            File circosFile = m.getCircos(cross, sample, accession);

            LineReader lr = new LineReader(circosFile);
            StringBuilder circosXml = new StringBuilder();
            String line;
            while ((line = lr.getNextRecord()) != null) {
                if (!line.contains("id=\"contig")) {
                    circosXml.append(line).append("\n");
                }
            }

            write(httpExchange, circosXml.toString());
        }
    }

    private class Search extends BaseHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = query(httpExchange.getRequestURI().getQuery());

            String cross = query.get("cross");
            String sample = query.get("sample");
            String accession = query.get("accession");
            String jexl = query.get("jexl");

            Expression e = je.createExpression(jexl);

            Set<String> contigNames = new HashSet<String>();

            File metrics = m.getMetrics(cross, sample, accession);
            TableReader tr = new TableReader(metrics);

            for (Map<String, String> te : tr) {
                JexlContext jexlContext = new MapContext();

                for (String key : te.keySet()) {
                    jexlContext.set(key, te.get(key));
                }

                try {
                    if (!e.getExpression().isEmpty() && (Boolean) e.evaluate(jexlContext)) {
                        String contigName = te.get("contigName");
                        contigNames.add(contigName);
                    }
                } catch (ClassCastException ex) {
                    throw new IndianaException("Problem evaluating JEXL expression for expression " + e + ": ", ex);
                }
            }

            StringBuilder xml = new StringBuilder();

            Pattern p = Pattern.compile("id=\"(.+)\"");

            File circosFile = m.getCircos(cross, sample, accession);
            LineReader lr = new LineReader(circosFile);
            String line;
            while ((line = lr.getNextRecord()) != null) {
                if (line.contains("id=\"contig")) {
                    Matcher m = p.matcher(line);

                    if (m.find()) {
                        String contigName = m.group(1);

                        if (contigNames.contains(contigName)) {
                            xml.append(line);
                        }
                    }
                } else {
                    xml.append(line);
                }
            }

            write(httpExchange, xml.toString());
        }
    }

    private class Info extends BaseHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = query(httpExchange.getRequestURI().getQuery());

            String cross = query.get("cross");
            String sample = query.get("sample");
            String accession = query.get("accession");
            String contigName = query.get("contigName");

            File metrics = m.getMetrics(cross, sample, accession);
            TableReader tr = new TableReader(metrics);

            for (Map<String, String> te : tr) {
                if (contigName.equals(te.get("contigName"))) {
                    te.remove("kmerOrigin");
                    te.remove("seq");

                    JSONObject jo = new JSONObject();
                    jo.put("info", te);

                    write(httpExchange, jo.toString());
                    break;
                }
            }

            write(httpExchange, 404, "not found");
        }
    }

    private class Graphs extends BaseHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = query(httpExchange.getRequestURI().getQuery());

            String cross = query.get("cross");
            String sample = query.get("sample");
            String accession = query.get("accession");

            Set<Map<String, String>> entries = new HashSet<Map<String, String>>();

            Map<String, String> allEntry = new HashMap<String, String>();
            allEntry.put("kmerSize", "all");
            allEntry.put("label", "all");
            entries.add(allEntry);

            for (String kmerSize : m.getGraphKmerSizes(cross, sample, accession)) {
                for (String graphLabel : m.getGraphLabels(cross, sample, accession, kmerSize)) {
                    Map<String, String> entry = new HashMap<String, String>();

                    entry.put("kmerSize", kmerSize);
                    entry.put("label", graphLabel);
                    entries.add(entry);
                }
            }

            JSONObject jo = new JSONObject();
            jo.put("entries", entries);

            write(httpExchange, jo.toString());
        }
    }

    private class Contig extends BaseHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = query(httpExchange.getRequestURI().getQuery());

            String cross = query.get("cross");
            String sample = query.get("sample");
            String accession = query.get("accession");
            String contigName = query.get("contigName");

            File metrics = m.getMetrics(cross, sample, accession);
            TableReader tr = new TableReader(metrics);

            for (Map<String, String> te : tr) {
                if (contigName.equals(te.get("contigName"))) {
                    String seq = te.get("seq");

                    Set<Map<String, String>> entries = new HashSet<Map<String, String>>();

                    for (int i = 0; i < seq.length(); i++) {
                        Map<String, String> entry = new HashMap<String, String>();
                        String base = String.valueOf(seq.charAt(i));
                        String id = String.format("%s_%d", base, i);
                        entry.put("base", base);
                        entry.put("name", id);
                        entry.put("class", base);

                        if (i > 0) {
                            String prevBase = String.valueOf(seq.charAt(i - 1));
                            String prevId = String.format("%s_%d", prevBase, i - 1);
                            entry.put("parent", prevId);
                        }

                        entries.add(entry);
                    }

                    for (String kmerSize : m.getGraphKmerSizes(cross, sample, accession)) {
                        for (String label: m.getGraphLabels(cross, sample, accession, kmerSize)) {
                            Map<String, File> graphMap = m.getGraphs(cross, sample, accession, kmerSize, label);

                            for (String name : graphMap.keySet()) {
                                CortexGraph cg = new CortexGraph(graphMap.get(name));
                                int ks = cg.getKmerSize();

                                for (int i = 0; i <= seq.length() - ks; i++) {
                                    String kmer = seq.substring(i, i + ks);

                                    Set<String> outKmers = CortexUtils.getNextKmers(cg, kmer);

                                    if (i <= seq.length() - ks - 1) {
                                        String nextKmer = seq.substring(i + 1, i + 1 + ks);

                                        outKmers.remove(nextKmer);

                                        for (String outKmer : outKmers) {
                                            Map<String, String> entry = new HashMap<String, String>();

                                            String outEdge = String.valueOf(outKmer.charAt(outKmer.length() - 1));
                                            String id = String.format("%s_%d", outEdge, i + ks);

                                            entry.put("base", outEdge);
                                            entry.put("name", id);
                                            entry.put("class", "noncontig");
                                            entries.add(entry);

                                            String prevBase = String.valueOf(kmer.charAt(kmer.length() - 1));
                                            String prevId = String.format("%s_%d", prevBase, i + ks - 1);
                                            entry.put("parent", prevId);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    JSONObject jo = new JSONObject();
                    jo.put("seq", new JSONArray(entries));

                    write(httpExchange, jo.toString());
                    break;
                }
            }

            write(httpExchange, 404, "not found");
        }
    }

    private class LinearContig extends BaseHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = query(httpExchange.getRequestURI().getQuery());

            String cross = query.get("cross");
            String sample = query.get("sample");
            String accession = query.get("accession");
            String contigName = query.get("contigName");
            String kmerSize = query.get("kmerSize");
            String graphLabel = query.get("graphLabel");

            File metrics = m.getMetrics(cross, sample, accession);
            TableReader tr = new TableReader(metrics);
            for (Map<String, String> te : tr) {
                if (te.get("contigName").equals(contigName)) {
                    String contig = te.get("seq");
                    String kmerOrigin = te.get("kmerOrigin");

                    Map<Integer, Set<Map<String, String>>> bases = new HashMap<Integer, Set<Map<String, String>>>();

                    for (int i = 0; i < contig.length(); i++) {
                        Map<String, String> baseMap = new HashMap<String, String>();
                        String base = String.valueOf(contig.charAt(i));
                        String id = String.format("%s_%d", base, i);

                        baseMap.put("base", base);
                        baseMap.put("id", id);
                        baseMap.put("class", base);
                        baseMap.put("pos", String.valueOf(i));
                        baseMap.put("type", "contig");
                        if (i < kmerOrigin.length()) {
                            baseMap.put("annotation", String.valueOf(kmerOrigin.charAt(i)));
                        }

                        bases.put(i, new HashSet<Map<String, String>>());
                        bases.get(i).add(baseMap);
                    }

                    for (String ks : m.getGraphKmerSizes(cross, sample, accession)) {
                        if (ks.equals(kmerSize) || kmerSize.equals("all")) {
                            for (String gl : m.getGraphLabels(cross, sample, accession, ks)) {
                                if (gl.equals(graphLabel) || graphLabel.equals("all")) {
                                    Map<String, File> graphMap = m.getGraphs(cross, sample, accession, ks, gl);

                                    for (String label : graphMap.keySet()) {
                                        CortexGraph cg = new CortexGraph(graphMap.get(label));

                                        for (int i = 0; i <= contig.length() - cg.getKmerSize(); i++) {
                                            String sk = contig.substring(i, i + cg.getKmerSize());
                                            CortexKmer ck = new CortexKmer(sk);
                                            CortexRecord cr = cg.findRecord(ck);

                                            if (cr == null) {
                                                for (Map<String, String> baseMap : bases.get(i)) {
                                                    baseMap.put("missing", "1");
                                                }
                                            } else {
                                                if (i > 0) {
                                                    String psk = contig.substring(i - 1, i - 1 + cg.getKmerSize());

                                                    Set<String> prevKmers = CortexUtils.getPrevKmers(cg, sk);
                                                    prevKmers.remove(psk);

                                                    if (prevKmers.size() > 0) {
                                                        for (String prevKmer : prevKmers) {
                                                            Map<String, String> prevBaseMap = new HashMap<String, String>();
                                                            String prevBase = String.valueOf(prevKmer.charAt(0));
                                                            String prevId = String.format("%s_%d", prevBase, i - 1);

                                                            prevBaseMap.put("base", prevBase);
                                                            prevBaseMap.put("id", prevId);
                                                            prevBaseMap.put("class", "none");
                                                            prevBaseMap.put("pos", String.valueOf(i - 1));
                                                            prevBaseMap.put("type", "inedge");

                                                            bases.get(i - 1).add(prevBaseMap);
                                                        }
                                                    }
                                                }

                                                if (i < contig.length() - cg.getKmerSize() - 1) {
                                                    String nsk = contig.substring(i + 1, i + 1 + cg.getKmerSize());

                                                    Set<String> nextKmers = CortexUtils.getNextKmers(cg, sk);
                                                    nextKmers.remove(nsk);

                                                    if (nextKmers.size() > 0) {
                                                        for (String nextKmer : nextKmers) {
                                                            Map<String, String> nextBaseMap = new HashMap<String, String>();
                                                            String nextBase = String.valueOf(nextKmer.charAt(nextKmer.length() - 1));
                                                            String nextId = String.format("%s_%d", nextBase, i + 1 + cg.getKmerSize());

                                                            nextBaseMap.put("base", nextBase);
                                                            nextBaseMap.put("id", nextId);
                                                            nextBaseMap.put("class", "none");
                                                            nextBaseMap.put("pos", String.valueOf(i + 1 + cg.getKmerSize()));
                                                            nextBaseMap.put("type", "outedge");

                                                            bases.get(i + 1 + cg.getKmerSize()).add(nextBaseMap);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    JSONObject jo = new JSONObject();
                    jo.put("bases", bases);
                    write(httpExchange, jo.toString());
                }
            }

            write(httpExchange, 404, "not implemented");
        }
    }

    @Override
    public void execute() {
        log.info("Loading...");

        m = new Manifest(MANIFEST);
        je = initializeJexlEngine();

        for (String cross : m.getCrosses()) {
            log.info("  {}: {} samples, {} accessions, {} kmer sizes, {} graphs, {} links",
                    cross,
                    m.numSamples(cross),
                    m.numAccessions(cross),
                    m.numKmerSizes(cross),
                    m.numGraphs(cross),
                    m.numLinks(cross)
            );
        }

        log.info("Starting server...");

        try {
            HttpServer server = HttpServer.create(new InetSocketAddress(PORT), 0);
            server.createContext("/", new PageHandler());
            server.createContext("/crosses", new CrossesList());
            server.createContext("/samples", new SamplesList());
            server.createContext("/accessions", new AccessionsList());
            server.createContext("/numcontigs", new NumContigs());
            server.createContext("/circos", new Circos());
            server.createContext("/search", new Search());
            server.createContext("/info", new Info());
            server.createContext("/graphs", new Graphs());
            server.createContext("/contig", new Contig());
            server.createContext("/linearcontig", new LinearContig());
            server.setExecutor(null);
            server.start();
        } catch (IOException e) {
            throw new IndianaException("Unable to start server", e);
        }

        log.info("  listening on port {}", PORT);

        /*
        contigs = loadContigs();
        links = loadLinks();

        log.info("  loaded {} contigs", contigs.size());
        for (String graphLabel : links.keySet()) {
            log.info("  loaded {} kmers with links from {}", links.get(graphLabel).size(), graphLabel);
        }

        log.info("Starting server...");

        try {
            HttpServer server = HttpServer.create(new InetSocketAddress(PORT), 0);
            server.createContext("/",                new PageHandler("/html/index.html"));
            server.createContext("/test",            new PageHandler("/html/newindex.html"));
            server.createContext("/tidy",            new PageHandler("/html/tidy.html"));
            server.createContext("/d3.v3.min.js",    new PageHandler("/html/d3.v3.min.js"));
            server.createContext("/autocomplete.js", new PageHandler("/html/autocomplete.js"));
            server.createContext("/indiana.css",     new PageHandler("/html/indiana.css"));
            server.createContext("/graph.json",      new PageHandler("/html/graph.json"));
            server.createContext("/contigs.csv",     new ContigsHandler());
            server.createContext("/contig",          new ContigHandler());
            server.createContext("/context",         new ContigContextHandler());
            server.createContext("/masked",          new MaskHandler());
            server.createContext("/edges",           new EdgeHandler());
            server.createContext("/links",           new LinksHandler());
            server.createContext("/link",            new LinkHandler());
            server.createContext("/linksplot",       new LinksPlotHandler());
            server.createContext("/covplot",         new CoveragePlotHandler());
            server.createContext("/dataset",         new DataSetHandler());
            server.createContext("/simplegraph",     new SimpleGraphHandler());
            server.setExecutor(null);
            server.start();
        } catch (IOException e) {
            throw new IndianaException("Unable to start server", e);
        }

        log.info("  listening on port {}", PORT);
        */
    }

    /*
    private Map<String, String> contigs;
    private Map<String, Map<CortexKmer, CortexLinksRecord>> links;

    private class PageHandler implements HttpHandler {
        private String page;

        public PageHandler(String page) {
            this.page = page;
        }

        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            InputStream is;

            if (log.isDebugEnabled()) {
                is = new FileInputStream(new File("./public" + page));

                log.debug("Reading '{}' from disk", page);
            } else {
                is = this.getClass().getResourceAsStream(page);
            }

            BufferedReader in = new BufferedReader(new InputStreamReader(is));
            StringBuilder responseData = new StringBuilder();

            String line;
            while((line = in.readLine()) != null) {
                responseData.append(line).append("\n");
            }

            String response = responseData.toString();
            response = response.replaceAll("\\$NUM_CONTIGS", String.valueOf(contigs.size()));

            if (page.endsWith(".js")) {
                Headers h = httpExchange.getResponseHeaders();
                h.add("Content-Type", "text/javascript");
            }

            httpExchange.sendResponseHeaders(200, response.length());

            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("Fetch static page : {}; response length: {}", httpExchange.getRequestURI(), response.length());
        }
    }

    private class ContigsHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            StringBuilder rb = new StringBuilder();
            rb.append("\"contigName\",\"baseLength\"\n");
            for (String name : contigs.keySet()) {
                rb.append("\"").append(name).append("\",\"").append(contigs.get(name).length()).append("\"\n");
            }

            String response = rb.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET contigs list  : {}; response length: {}", httpExchange.getRequestURI(), response.length());
        }
    }

    private Map<String, String> queryToMap(String query){
        Map<String, String> result = new HashMap<String, String>();
        for (String param : query.split("&")) {
            String pair[] = param.split("=");
            if (pair.length>1) {
                result.put(pair[0], pair[1]);
            }else{
                result.put(pair[0], "");
            }
        }
        return result;
    }

    private class ContigHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = queryToMap(httpExchange.getRequestURI().getQuery());

            StringBuilder rb = new StringBuilder();
            rb.append("{");
            rb.append("\"seq\": \"").append(contigs.get(query.get("contigName"))).append("\"");
            rb.append("}\n");

            String response = rb.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET contig seq    : {}, response length: {}", query.get("contigName"), response.length());
        }
    }

    private class ContigContextHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = queryToMap(httpExchange.getRequestURI().getQuery());

            String seq = contigs.get(query.get("contigName"));
            String selectedGraph = query.get("graphName");

            Map<Integer, Set<String>> contigContext = new TreeMap<Integer, Set<String>>();

            int earliestPosition = 0;
            int latestPosition = 0;

            for (String graphLabel : GRAPHS.keySet()) {
                if (selectedGraph.equals("all") || graphLabel.equals(selectedGraph)) {
                    CortexGraph cg = new CortexGraph(GRAPHS.get(graphLabel));
                    int kmerSize = cg.getKmerSize();

                    String currentKmer = seq.substring(0, kmerSize);
                    int currentPos = 0;

                    while (currentKmer != null) {
                        Set<String> prevKmers = CortexUtils.getPrevKmers(cg, currentKmer);

                        if (prevKmers.size() == 1) {
                            String prevKmer = prevKmers.iterator().next();

                            if (!contigContext.containsKey(currentPos - 1)) {
                                contigContext.put(currentPos - 1, new TreeSet<String>());
                            }

                            contigContext.get(currentPos - 1).add(String.valueOf(prevKmer.charAt(0)));

                            if (currentPos - 1 < earliestPosition) { earliestPosition = currentPos - 1; }

                            currentKmer = prevKmer;
                            currentPos--;
                        } else {
                            currentKmer = null;
                        }
                    }

                    currentKmer = seq.substring(seq.length() - kmerSize, seq.length());
                    currentPos = seq.length() - 1;

                    while (currentKmer != null) {
                        Set<String> nextKmers = CortexUtils.getNextKmers(cg, currentKmer);

                        if (nextKmers.size() == 1) {
                            String nextKmer = nextKmers.iterator().next();

                            if (!contigContext.containsKey(currentPos + 1)) {
                                contigContext.put(currentPos + 1, new TreeSet<String>());
                            }

                            contigContext.get(currentPos + 1).add(String.valueOf(nextKmer.charAt(nextKmer.length() - 1)));

                            if (currentPos + 1 > latestPosition) { latestPosition = currentPos + 1; }

                            currentKmer = nextKmer;
                            currentPos++;
                        } else {
                            currentKmer = null;
                        }
                    }
                }
            }

            for (int i = -1; i >= earliestPosition; i--) {
                if (contigContext.containsKey(i) && contigContext.get(i).size() > 1) {
                    earliestPosition = i + 1;
                    break;
                }
            }

            for (int i = seq.length(); i <= latestPosition; i++) {
                if (contigContext.containsKey(i) && contigContext.get(i).size() > 1) {
                    latestPosition = i - 1;
                    break;
                }
            }

            Map<Integer, String> contextMap = new TreeMap<Integer, String>();

            for (Integer pos : contigContext.keySet()) {
                if (pos >= earliestPosition && pos <= latestPosition) {
                    contextMap.put(pos, contigContext.get(pos).iterator().next());
                }
            }

            JSONObject jo = new JSONObject(contextMap);

            String response = jo.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET context       : {}, response length: {}", query.get("contigName"), response.length());
        }
    }

    private class JSONTree {
        private String root;
        private Set<JSONTree> children;

        public JSONTree(char root) {
            this.root = String.valueOf(root);
        }

        public JSONTree(String root) {
            this.root = root;
        }

        public JSONTree add(char child) {
            return add(String.valueOf(child));
        }

        public JSONTree add(String child) {
            if (children == null) {
                children = new HashSet<JSONTree>();
            }

            JSONTree childTree = new JSONTree(child);

            children.add(childTree);

            return childTree;
        }

        public String getRoot() { return root; }

        private StringBuilder traverse() {
            Collection<StringBuilder> co = new ArrayList<StringBuilder>();
            StringBuilder jo = new StringBuilder();
            //jo.put("base", root);
            jo.append("{\"base\": \"").append(root).append("\"");

            if (children != null) {
                List<StringBuilder> pieces = new ArrayList<StringBuilder>();
                for (JSONTree childTree : children) {
                    pieces.add(childTree.traverse());
                }

                jo.append(", \"children\": [").append(Joiner.on(",").join(pieces)).append("]");

                //jo.put("children", new JSONArray(co));
                jo.append("},");
            } else {
                jo.append("}");
            }

            return jo;
        }

        public String toString() {
            StringBuilder jo = traverse();

            return jo.toString();
        }
    }

    private class SimpleGraphHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = queryToMap(httpExchange.getRequestURI().getQuery());

            //String seq = contigs.get(query.get("contigName"));
            String seq = "TACG";
            String selectedGraph = query.get("graphName");

            JSONTree jtree = null;

            for (String graphLabel : GRAPHS.keySet()) {
                if (selectedGraph.equals("all") || graphLabel.equals(selectedGraph)) {
                    //CortexGraph cg = new CortexGraph(GRAPHS.get(graphLabel));
                    //int kmerSize = cg.getKmerSize();
                    int kmerSize = 2;

                    JSONTree jroot = null;

                    for (int i = 0; i <= seq.length() - kmerSize; i++) {
                        String kmer = seq.substring(i, i + kmerSize);

                        for (int j = 0; j < kmer.length(); j++) {
                            if (jtree == null) {
                                jtree = new JSONTree(kmer.charAt(0));
                                jroot = jtree;
                            } else {
                                if (jroot != null) {
                                    jroot = jroot.add(kmer.charAt(j));
                                }
                            }
                        }
                    }

                    log.info("jtree: {}", jtree);
                }
            }

            String response = jtree == null ? "empty" : jtree.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            //log.info("GET context       : {}, response length: {}", query.get("contigName"), response.length());
        }
    }

    private class EdgeHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = queryToMap(httpExchange.getRequestURI().getQuery());

            String seq = contigs.get(query.get("contigName"));
            String selectedGraph = query.get("graphName");

            Map<Integer, Map<String, Set<String>>> allEdges = new HashMap<Integer, Map<String, Set<String>>>();

            for (String graphLabel : GRAPHS.keySet()) {
                if (selectedGraph.equals("all") || graphLabel.equals(selectedGraph)) {
                    CortexGraph cg = new CortexGraph(GRAPHS.get(graphLabel));
                    int kmerSize = cg.getKmerSize();

                    String currentKmer = seq.substring(0, kmerSize);
                    int currentPos = 0;

                    while (currentKmer != null) {
                        Set<String> prevKmers = CortexUtils.getPrevKmers(cg, currentKmer);

                        if (prevKmers.size() == 1) {
                            currentKmer = prevKmers.iterator().next();
                            currentPos--;
                        } else {
                            int iePos = currentPos - 1;

                            if (prevKmers.size() > 1) {
                                if (!allEdges.containsKey(iePos)) {
                                    allEdges.put(iePos, new HashMap<String, Set<String>>());
                                }

                                if (!allEdges.get(iePos).containsKey("in")) {
                                    allEdges.get(iePos).put("in", new TreeSet<String>());
                                }

                                for (String prevKmer : prevKmers) {
                                    allEdges.get(iePos).get("in").add(prevKmer.substring(0, 1));
                                }
                            }

                            currentKmer = null;
                        }
                    }

                    currentKmer = seq.substring(seq.length() - kmerSize, seq.length());
                    currentPos = seq.length() - 1;

                    while (currentKmer != null) {
                        Set<String> nextKmers = CortexUtils.getNextKmers(cg, currentKmer);

                        if (nextKmers.size() == 1) {
                            currentKmer = nextKmers.iterator().next();
                            currentPos++;
                        } else {
                            int oePos = currentPos - 1;

                            if (nextKmers.size() > 1) {
                                if (!allEdges.containsKey(oePos)) {
                                    allEdges.put(oePos, new HashMap<String, Set<String>>());
                                }

                                if (!allEdges.get(oePos).containsKey("out")) {
                                    allEdges.get(oePos).put("out", new TreeSet<String>());
                                }

                                for (String nextKmer : nextKmers) {
                                    allEdges.get(oePos).get("out").add(nextKmer.substring(nextKmer.length() - 1, nextKmer.length()));
                                }
                            }
                            currentKmer = null;
                        }
                    }

                    for (int i = 1; i <= seq.length() - kmerSize - 1; i++) {
                        int iePos = i - 1;
                        int oePos = i + kmerSize;

                        String kmer = seq.substring(i, i + kmerSize);
                        CortexKmer ck = new CortexKmer(kmer);
                        CortexRecord cr = cg.findRecord(ck);

                        if (cr != null) {
                            Set<String> ie = CortexUtils.getPrevKmers(cg, kmer);
                            Set<String> oe = CortexUtils.getNextKmers(cg, kmer);

                            if (i > 0) {
                                String prevKmer = seq.substring(iePos, oePos - 1);
                                ie.remove(prevKmer);
                            }

                            if (i < seq.length() - kmerSize) {
                                String nextKmer = seq.substring(i + 1, i + 1 + kmerSize);
                                oe.remove(nextKmer);
                            }

                            if (ie.size() > 0) {
                                if (!allEdges.containsKey(iePos)) {
                                    allEdges.put(iePos, new HashMap<String, Set<String>>());
                                }

                                if (!allEdges.get(iePos).containsKey("in")) {
                                    allEdges.get(iePos).put("in", new TreeSet<String>());
                                }

                                for (String inKmer : ie) {
                                    allEdges.get(iePos).get("in").add(inKmer.substring(0, 1));
                                }
                            }

                            if (oe.size() > 0) {
                                if (!allEdges.containsKey(oePos)) {
                                    allEdges.put(oePos, new HashMap<String, Set<String>>());
                                }

                                if (!allEdges.get(oePos).containsKey("out")) {
                                    allEdges.get(oePos).put("out", new TreeSet<String>());
                                }

                                for (String outKmer : oe) {
                                    allEdges.get(oePos).get("out").add(outKmer.substring(outKmer.length() - 1, outKmer.length()));
                                }
                            }
                        }
                    }
                }
            }

            String response = (new JSONObject(allEdges)).toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET in/out edges  : {}, response length: {}", query.get("contigName"), response.length());
        }
    }

    private class MaskHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = queryToMap(httpExchange.getRequestURI().getQuery());

            String seq = contigs.get(query.get("contigName"));
            String selectedGraph = query.get("graphName");

            Set<Integer> maskedPositions = new TreeSet<Integer>();

            for (String graphLabel : GRAPHS.keySet()) {
                if (selectedGraph.equals("all") || graphLabel.equals(selectedGraph)) {
                    CortexGraph cg = new CortexGraph(GRAPHS.get(graphLabel));
                    int kmerSize = cg.getKmerSize();

                    for (int i = 0; i <= seq.length() - kmerSize; i++) {
                        String thisKmer = seq.substring(i, i + kmerSize);
                        CortexKmer thisCk = new CortexKmer(thisKmer);
                        CortexRecord thisCr = cg.findRecord(thisCk);

                        if (thisCr == null) {
                            maskedPositions.add(i);
                        }
                    }
                }
            }

            JSONObject jsonResponse = new JSONObject();
            jsonResponse.put("masked", new JSONArray(maskedPositions));

            String response = jsonResponse.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET masked bases  : current, response length: {}", response.length());
        }
    }

    private class LinksHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = queryToMap(httpExchange.getRequestURI().getQuery());

            String seq = contigs.get(query.get("contigName"));
            String selectedGraph = query.get("graphName");
            boolean goForward = query.get("orientation").equals("forward");

            Set<Integer> linkStarts = new TreeSet<Integer>();

            for (String graphLabel : LINKS.keySet()) {
                if (selectedGraph.equals("all") || graphLabel.equals(selectedGraph)) {
                    Map<CortexKmer, CortexLinksRecord> cprs = links.get(graphLabel);
                    CortexGraph cg = new CortexGraph(GRAPHS.get(graphLabel));
                    int kmerSize = cg.getKmerSize();

                    for (int i = 0; i <= seq.length() - kmerSize; i++) {
                        String thisKmer = seq.substring(i, i + kmerSize);
                        CortexKmer thisCk = new CortexKmer(thisKmer);

                        if (cprs.containsKey(thisCk)) {
                            CortexLinksRecord clr = cprs.get(thisCk);

                            for (CortexJunctionsRecord cjr : clr.getJunctions()) {
                                if (thisCk.isFlipped()) {
                                    cjr = CortexUtils.flipJunctionsRecord(cjr);
                                }

                                if (goForward == cjr.isForward()) {
                                    if (goForward) {
                                        linkStarts.add(i);
                                    } else {
                                        linkStarts.add(i + kmerSize - 1);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            JSONObject jsonResponse = new JSONObject();
            jsonResponse.put("linkStarts", new JSONArray(linkStarts));

            String response = jsonResponse.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET links info    : current, response length: {}", response.length());
        }
    }

    private class LinkEntry {
        public String linkStartId;
        public String linkEndId;

        public LinkEntry(String linkStartId, String linkEndId) {
            this.linkStartId = linkStartId;
            this.linkEndId = linkEndId;
        }

        public Map<String, String> getJSONMap() {
            Map<String, String> jsonMap = new LinkedHashMap<String, String>();
            jsonMap.put("linkStartId", linkStartId);
            jsonMap.put("linkEndId", linkEndId);

            return jsonMap;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            LinkEntry linkEntry = (LinkEntry) o;

            if (linkEndId != null ? !linkEndId.equals(linkEntry.linkEndId) : linkEntry.linkEndId != null) return false;
            if (linkStartId != null ? !linkStartId.equals(linkEntry.linkStartId) : linkEntry.linkStartId != null)
                return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = linkStartId != null ? linkStartId.hashCode() : 0;
            result = 31 * result + (linkEndId != null ? linkEndId.hashCode() : 0);
            return result;
        }
    }

    private class CoveragePlotHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = queryToMap(httpExchange.getRequestURI().getQuery());

            String seq = contigs.get(query.get("contigName"));
            String selectedGraph = query.get("graphName");

            Map<Integer, Integer> kmerCoverage = new TreeMap<Integer, Integer>();

            for (int pos = 0; pos <= seq.length(); pos++) {
                kmerCoverage.put(pos, 0);
            }

            for (String graphLabel : GRAPHS.keySet()) {
                if (selectedGraph.equals("all") || graphLabel.equals(selectedGraph)) {
                    CortexGraph cg = new CortexGraph(GRAPHS.get(graphLabel));
                    int kmerSize = cg.getKmerSize();

                    for (int pos = 0; pos <= seq.length() - kmerSize; pos++) {
                        String sk = seq.substring(pos, pos + kmerSize);
                        CortexKmer ck = new CortexKmer(sk);
                        CortexRecord cr = cg.findRecord(ck);

                        if (cr != null) {
                            kmerCoverage.put(pos, kmerCoverage.get(pos) + cr.getCoverage(0));
                        }
                    }
                }
            }

            JSONObject jsonResponse = new JSONObject();
            jsonResponse.put("cov", kmerCoverage.values());

            String response = jsonResponse.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET kmer cov info : {}, response length: {}", query.get("contigName"), response.length());
        }
    }

    private class LinksPlotHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = queryToMap(httpExchange.getRequestURI().getQuery());

            String seq = contigs.get(query.get("contigName"));
            String selectedGraph = query.get("graphName");

            Map<Integer, Integer> linkCoverage = new TreeMap<Integer, Integer>();

            for (String graphLabel : LINKS.keySet()) {
                if (selectedGraph.equals("all") || graphLabel.equals(selectedGraph)) {
                    Map<CortexKmer, CortexLinksRecord> clrs = links.get(graphLabel);
                    CortexGraph cg = new CortexGraph(GRAPHS.get(graphLabel));
                    int kmerSize = cg.getKmerSize();

                    for (int pos = 0; pos <= seq.length(); pos++) {
                        linkCoverage.put(pos, 0);
                    }

                    for (int pos = 0; pos <= seq.length() - kmerSize; pos++) {
                        String sk = seq.substring(pos, pos + kmerSize);
                        CortexKmer ck = new CortexKmer(sk);

                        if (!clrs.containsKey(ck) && pos - kmerSize + 1 >= 0) {
                            sk = seq.substring(pos - kmerSize + 1, pos + 1);
                            ck = new CortexKmer(sk);
                        }

                        if (clrs.containsKey(ck)) {
                            for (CortexJunctionsRecord cjr : clrs.get(ck).getJunctions()) {
                                if (ck.isFlipped()) {
                                    cjr = CortexUtils.flipJunctionsRecord(cjr);
                                }

                                String junctions = cjr.getJunctions();

                                if (cjr.isForward()) {
                                    int currentJunction = 0;
                                    int startPos = pos;

                                    for (int i = pos; i <= seq.length() - kmerSize && currentJunction < junctions.length(); i++) {
                                        String thisKmer = seq.substring(i, i + kmerSize);
                                        Set<String> nextKmers = CortexUtils.getNextKmers(cg, thisKmer);
                                        if (nextKmers.size() > 1) {
                                            String expectedNextKmer = seq.substring(i + 1, i + kmerSize) + junctions.charAt(currentJunction);

                                            if (nextKmers.contains(expectedNextKmer)) {
                                                int endPos = i + kmerSize;

                                                for (int j = startPos; j <= endPos; j++) {
                                                    linkCoverage.put(j, linkCoverage.get(j) + 1);
                                                }

                                                startPos = endPos;

                                                currentJunction++;

                                                if (i < seq.length() - kmerSize) {
                                                    String nextKmerInContig = seq.substring(i + 1, i + 1 + kmerSize);

                                                    if (!expectedNextKmer.equals(nextKmerInContig)) {
                                                        break;
                                                    }
                                                }
                                            } else {
                                                //log.info("Did not find expected next kmer: {} {}", i, expectedNextKmer);
                                            }
                                        } else if (nextKmers.size() == 0) {
                                            break;
                                        }
                                    }
                                } else {
                                    int currentJunction = 0;
                                    int startPos = pos;

                                    for (int i = pos; i >= kmerSize && currentJunction < junctions.length(); i--) {
                                        String thisKmer = seq.substring(i - kmerSize + 1, i + 1);
                                        Set<String> prevKmers = CortexUtils.getPrevKmers(cg, thisKmer);
                                        if (prevKmers.size() > 1) {
                                            String expectedPrevKmer = junctions.charAt(currentJunction) + seq.substring(i - kmerSize + 1, i);

                                            if (prevKmers.contains(expectedPrevKmer)) {
                                                int endPos = i - kmerSize;

                                                for (int j = startPos; j <= endPos; j++) {
                                                    linkCoverage.put(j, linkCoverage.get(j) + 1);
                                                }

                                                startPos = endPos;

                                                currentJunction++;

                                                if (i > kmerSize) {
                                                    String prevKmerInContig = seq.substring(i - kmerSize, i);

                                                    if (!expectedPrevKmer.equals(prevKmerInContig)) {
                                                        break;
                                                    }
                                                }
                                            } else {
                                                //log.info("Did not find expected prev kmer: {} {}", i, expectedPrevKmer);
                                            }
                                        } else if (prevKmers.size() == 0) {
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            JSONObject jsonResponse = new JSONObject();
            jsonResponse.put("cov", linkCoverage.values());

            String response = jsonResponse.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET link plot info: current, response length: {}", response.length());
        }
    }

    private class LinkHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            Map<String, String> query = queryToMap(httpExchange.getRequestURI().getQuery());

            String seq = contigs.get(query.get("contigName"));
            String selectedGraph = query.get("graphName");
            int pos = Integer.valueOf(query.get("pos"));

            Map<LinkEntry, Integer> linkEntries = new LinkedHashMap<LinkEntry, Integer>();

            for (String graphLabel : LINKS.keySet()) {
                if (selectedGraph.equals("all") || graphLabel.equals(selectedGraph)) {
                    Map<CortexKmer, CortexLinksRecord> clrs = links.get(graphLabel);
                    CortexGraph cg = new CortexGraph(GRAPHS.get(graphLabel));
                    int kmerSize = cg.getKmerSize();

                    String sk = seq.substring(pos, pos + kmerSize);
                    CortexKmer ck = new CortexKmer(sk);

                    if (!clrs.containsKey(ck)) {
                        sk = seq.substring(pos - kmerSize + 1, pos + 1);
                        ck = new CortexKmer(sk);
                    }

                    for (CortexJunctionsRecord cjr : clrs.get(ck).getJunctions()) {
                        if (ck.isFlipped()) {
                            cjr = CortexUtils.flipJunctionsRecord(cjr);
                        }

                        String junctions = cjr.getJunctions();

                        if (cjr.isForward()) {
                            log.info("junctions: {}", cjr);

                            int currentJunction = 0;
                            int startPos = pos;
                            String startBase = sk.substring(0, 1);

                            for (int i = pos; i <= seq.length() - kmerSize && currentJunction < junctions.length(); i++) {
                                String thisKmer = seq.substring(i, i + kmerSize);
                                Set<String> nextKmers = CortexUtils.getNextKmers(cg, thisKmer);
                                if (nextKmers.size() > 1) {
                                    String expectedNextKmer = seq.substring(i + 1, i + kmerSize) + junctions.charAt(currentJunction);

                                    if (nextKmers.contains(expectedNextKmer)) {
                                        int endPos = i + kmerSize;
                                        String endBase = junctions.substring(currentJunction, currentJunction + 1);

                                        LinkEntry le1 = new LinkEntry(
                                            String.format("%s_%d", startBase, startPos),
                                            String.format("%s_%d", expectedNextKmer.substring(expectedNextKmer.length() - 2, expectedNextKmer.length() - 1), i + kmerSize - 1));

                                        LinkEntry le2 = new LinkEntry(
                                            String.format("%s_%d", expectedNextKmer.substring(expectedNextKmer.length() - 2, expectedNextKmer.length() - 1), i + kmerSize - 1),
                                            String.format("%s_%d", endBase, endPos));

                                        if (linkEntries.containsKey(le1)) {
                                            linkEntries.put(le1, linkEntries.get(le1) + cjr.getCoverage(0));
                                        } else {
                                            linkEntries.put(le1, cjr.getCoverage(0));
                                        }

                                        if (linkEntries.containsKey(le2)) {
                                            linkEntries.put(le2, linkEntries.get(le2) + cjr.getCoverage(0));
                                        } else {
                                            linkEntries.put(le2, cjr.getCoverage(0));
                                        }

                                        startPos = endPos;
                                        startBase = endBase;

                                        currentJunction++;

                                        if (i < seq.length() - kmerSize) {
                                            String nextKmerInContig = seq.substring(i + 1, i + 1 + kmerSize);

                                            if (!expectedNextKmer.equals(nextKmerInContig)) {
                                                break;
                                            }
                                        }
                                    } else {
                                        log.info("Did not find expected next kmer: {} {}", i, expectedNextKmer);
                                    }
                                } else if (nextKmers.size() == 0) {
                                    log.info("Could not navigate forward from kmer: {} {}", i, thisKmer);
                                    break;
                                }
                            }
                        } else {
                            log.info("junctions: {}", cjr);

                            int currentJunction = 0;
                            int startPos = pos;
                            String startBase = sk.substring(sk.length() - 1);

                            for (int i = pos; i >= kmerSize && currentJunction < junctions.length(); i--) {
                                String thisKmer = seq.substring(i - kmerSize + 1, i + 1);
                                Set<String> prevKmers = CortexUtils.getPrevKmers(cg, thisKmer);
                                if (prevKmers.size() > 1) {
                                    String expectedPrevKmer = junctions.charAt(currentJunction) + seq.substring(i - kmerSize + 1, i);

                                    if (prevKmers.contains(expectedPrevKmer)) {
                                        int endPos = i - kmerSize;
                                        String endBase = junctions.substring(currentJunction, currentJunction + 1);

                                        LinkEntry le1 = new LinkEntry(
                                                String.format("%s_%d", startBase, startPos),
                                                String.format("%s_%d", expectedPrevKmer.substring(1, 2), i - kmerSize + 1));

                                        LinkEntry le2 = new LinkEntry(
                                                String.format("%s_%d", expectedPrevKmer.substring(1, 2), i - kmerSize + 1),
                                                String.format("%s_%d", endBase, endPos));

                                        if (linkEntries.containsKey(le1)) {
                                            linkEntries.put(le1, linkEntries.get(le1) + cjr.getCoverage(0));
                                        } else {
                                            linkEntries.put(le1, cjr.getCoverage(0));
                                        }

                                        if (linkEntries.containsKey(le2)) {
                                            linkEntries.put(le2, linkEntries.get(le2) + cjr.getCoverage(0));
                                        } else {
                                            linkEntries.put(le2, cjr.getCoverage(0));
                                        }

                                        startPos = endPos;
                                        startBase = endBase;

                                        currentJunction++;

                                        if (i > kmerSize) {
                                            String prevKmerInContig = seq.substring(i - kmerSize, i);

                                            if (!expectedPrevKmer.equals(prevKmerInContig)) {
                                                break;
                                            }
                                        }
                                    } else {
                                        log.info("Did not find expected prev kmer: {} {}", i, expectedPrevKmer);
                                    }
                                } else if (prevKmers.size() == 0) {
                                    log.info("Could not navigate backward from kmer: {} {}", i, thisKmer);
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            JSONObject jsonResponse = new JSONObject();

            for (LinkEntry le : linkEntries.keySet()) {
                Map<String, String> jsonMap = le.getJSONMap();
                jsonMap.put("linkCoverage", String.valueOf(linkEntries.get(le)));

                jsonResponse.append("linkEntries", jsonMap);
            }

            String response = jsonResponse.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET link info     : current, response length: {}", response.length());
        }
    }

    private class DataSetHandler implements HttpHandler {
        @Override
        public void handle(HttpExchange httpExchange) throws IOException {
            JSONObject jsonResponse = new JSONObject();
            jsonResponse.put("contigsFile", CONTIGS.getAbsolutePath());

            JSONArray graphInfos = new JSONArray();
            for (String graphName : GRAPHS.keySet()) {
                Map<String, Object> graphInfo = new HashMap<String, Object>();
                CortexGraph cg = new CortexGraph(GRAPHS.get(graphName));

                graphInfo.put("graphName", graphName);
                graphInfo.put("graphFile", GRAPHS.get(graphName).getAbsolutePath());
                graphInfo.put("kmerSize", cg.getKmerSize());
                graphInfo.put("numRecords", cg.getNumRecords());
                graphInfo.put("sampleName", cg.getColor(0).getSampleName());

                graphInfos.put(graphInfo);
            }

            jsonResponse.put("graphs", graphInfos);

            String response = jsonResponse.toString();

            httpExchange.sendResponseHeaders(200, response.length());
            OutputStream os = httpExchange.getResponseBody();
            os.write(response.getBytes());
            os.close();

            log.info("GET dataset info  : current, response length: {}", response.length());
        }
    }

    private Map<String, String> loadContigs() {
        FastaSequenceFile contigReader = new FastaSequenceFile(CONTIGS, true);

        Map<String, String> contigs = new LinkedHashMap<String, String>();

        ReferenceSequence rseq;
        while ((rseq = contigReader.nextSequence()) != null) {
            String name = rseq.getName();

            contigs.put(name, new String(rseq.getBases()));
        }

        return contigs;
    }

    private Map<String, Map<CortexKmer, CortexLinksRecord>> loadLinks() {
        Map<String, Map<CortexKmer, CortexLinksRecord>> links = new HashMap<String, Map<CortexKmer, CortexLinksRecord>>();

        for (String graphName : LINKS.keySet()) {
            CortexLinks cp = new CortexLinks(LINKS.get(graphName));

            links.put(graphName, new HashMap<CortexKmer, CortexLinksRecord>());

            for (CortexLinksRecord cpr : cp) {
                links.get(graphName).put(cpr.getBinaryKmer(), cpr);
            }
        }

        return links;
    }
    */

}
