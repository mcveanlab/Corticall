package uk.ac.ox.well.indiana.attic.analyses.nahr;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.StringUtil;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksMap;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.PrintStream;
import java.util.*;

public class ContigConfidence extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (.fasta)")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="graph", shortName="g", doc="Graph (.ctx)")
    public LinkedHashSet<CortexGraph> GRAPHS;

    @Argument(fullName="links", shortName="l", doc="Links (.ctp)")
    public HashSet<CortexLinksMap> LINKS;

    @Argument(fullName="contigName", shortName="n", doc="Contig name")
    public HashSet<String> CONTIG_NAMES;

    @Output
    public PrintStream out;

    private class GraphAndLinks {
        private CortexGraph cg;
        private CortexLinksMap clm;

        public void setGraph(CortexGraph cg) { this.cg = cg; }
        public void setLinks(CortexLinksMap clm) { this.clm = clm; }

        public CortexGraph getGraph() { return cg; }
        public CortexLinksMap getLinks() { return clm; }
    }

    @Override
    public void execute() {
        log.info("Loading contigs...");
        Map<String, ReferenceSequence> contigs = new HashMap<String, ReferenceSequence>();
        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            contigs.put(rseq.getName(), rseq);
        }
        log.info("  {} contigs", contigs.size());

        log.info("Loading graphs and links...");
        Map<String, GraphAndLinks> graphAndLinksMap = new HashMap<String, GraphAndLinks>();
        for (CortexGraph cg : GRAPHS) {
            String sample = cg.getColor(0).getSampleName() + ".k" + cg.getKmerSize();

            GraphAndLinks gl = new GraphAndLinks();
            gl.setGraph(cg);

            graphAndLinksMap.put(sample, gl);
        }

        for (CortexLinksMap clm : LINKS) {
            String sample = clm.getCortexLinks().getColor(0).getSampleName() + ".k" + clm.getCortexLinks().getKmerSize();

            if (graphAndLinksMap.containsKey(sample)) {
                graphAndLinksMap.get(sample).setLinks(clm);
            }
        }

        for (String sample : graphAndLinksMap.keySet()) {
            GraphAndLinks gl = graphAndLinksMap.get(sample);
            String cgName = (gl.getGraph() == null) ? "none" : gl.getGraph().getCortexFile().getName();
            String clName = (gl.getLinks() == null) ? "none" : gl.getLinks().getCortexLinks().getCortexLinksFile().getName();

            log.info("  {}: {} {}", sample, cgName, clName);
        }

        for (String contigName : CONTIG_NAMES) {
            log.info("contigName: {}", contigName);

            String seq = new String(contigs.get(contigName).getBases());

            for (String sample : graphAndLinksMap.keySet()) {
                GraphAndLinks gl = graphAndLinksMap.get(sample);

                if (gl.getGraph() != null && gl.getLinks() != null) {
                    CortexGraph cg = gl.getGraph();
                    CortexLinksMap clm = gl.getLinks();

                    log.info("Graph: {}", cg.getCortexFile().getName());

                    StringBuilder firstFlank = new StringBuilder();
                    String firstKmer = seq.substring(0, cg.getKmerSize());
                    Set<String> pks = CortexUtils.getPrevKmers(cg, firstKmer);
                    while (pks.size() == 1) {
                        String kmer = pks.iterator().next();
                        firstFlank.insert(0, kmer.charAt(0));

                        pks = CortexUtils.getPrevKmers(cg, kmer);
                    }

                    StringBuilder lastFlank = new StringBuilder();
                    String lastKmer = seq.substring(seq.length() - cg.getKmerSize(), seq.length());
                    Set<String> nks = CortexUtils.getNextKmers(cg, lastKmer);
                    while (nks.size() == 1) {
                        String kmer = nks.iterator().next();
                        lastFlank.append(kmer.charAt(kmer.length() - 1));

                        nks = CortexUtils.getNextKmers(cg, kmer);
                    }

                    String contig = firstFlank.toString() + seq + lastFlank.toString();

                    byte[] oedges = new byte[contig.length()];
                    byte[] iedges = new byte[contig.length()];

                    byte[][] oeb = new byte[3][contig.length()];
                    byte[][] ieb = new byte[3][contig.length()];

                    byte[] olink = new byte[contig.length()];
                    byte[] ilink = new byte[contig.length()];
                    byte[] xaxis = new byte[contig.length()];

                    for (int i = 0; i < contig.length(); i++) {
                        oedges[i] = ' ';
                        iedges[i] = ' ';

                        for (int j = 0; j < 3; j++) {
                            oeb[j][i] = ' ';
                            ieb[j][i] = ' ';
                        }

                        olink[i] = ' ';
                        ilink[i] = ' ';
                        xaxis[i] = (byte) ((i % 25 == 0) ? '.' : ' ');
                    }

                    Map<Integer, Set<String>> kmers = new TreeMap<Integer, Set<String>>();
                    Set<String> contigKmers = new HashSet<String>();

                    Map<String, CortexLinksRecord> linkMap = new HashMap<String, CortexLinksRecord>();
                    Map<CortexJunctionsRecord, Set<String>> kmersInLinks = new HashMap<CortexJunctionsRecord, Set<String>>();

                    for (int i = 0; i <= contig.length() - cg.getKmerSize(); i++) {
                        log.info("Base {}", i);

                        String sk = contig.substring(i, i + cg.getKmerSize());
                        CortexKmer ck = new CortexKmer(sk);

                        if (clm.containsKey(ck)) {
                            linkMap.put(sk, clm.get(ck));

                            for (CortexJunctionsRecord cjr : clm.get(ck).getJunctions()) {
                                Set<String> kil = new HashSet<String>(CortexUtils.getKmersInLinkByNavigation(cg, sk, cjr));

                                kmersInLinks.put(cjr, kil);
                            }
                        }
                    }

                    Set<String> iogt2 = new HashSet<String>();

                    for (int i = 0; i <= contig.length() - cg.getKmerSize(); i++) {
                        String sk = contig.substring(i, i + cg.getKmerSize());

                        Set<String> inKmers  = CortexUtils.getPrevKmers(cg, sk);
                        Set<String> outKmers = CortexUtils.getNextKmers(cg, sk);

                        if (inKmers.size() > 1) {
                            iedges[i] = '|';

                            char base = i > 0 ? contig.charAt(i - 1) : ' ';
                            for (String inKmer : inKmers) {
                                iogt2.add(inKmer);

                                if (inKmer.charAt(0) != base) {
                                    for (int j = 0; j < 3; j++) {
                                        if (ieb[j][i] == ' ') {
                                            ieb[j][i] = (byte) inKmer.charAt(0);
                                            break;
                                        }
                                    }
                                }
                            }
                        }

                        if (outKmers.size() > 1) {
                            oedges[i+cg.getKmerSize()-1] = '|';

                            char base = (i + cg.getKmerSize() < contig.length()) ? contig.charAt(i + cg.getKmerSize()) : ' ';
                            for (String outKmer : outKmers) {
                                iogt2.add(outKmer);

                                if (outKmer.charAt(outKmer.length() - 1) != base) {
                                    for (int j = 0; j < 3; j++) {
                                        if (oeb[j][i + cg.getKmerSize() - 1] == ' ') {
                                            oeb[j][i + cg.getKmerSize() - 1] = (byte) outKmer.charAt(outKmer.length() - 1);
                                            break;
                                        }
                                    }
                                }
                            }
                        }

                        if (!kmers.containsKey(i-1)) { kmers.put(i-1, new HashSet<String>()); }
                        if (!kmers.containsKey(i))   { kmers.put(i,   new HashSet<String>()); }
                        if (!kmers.containsKey(i+1)) { kmers.put(i+1, new HashSet<String>()); }

                        kmers.get(i - 1).addAll(inKmers);
                        kmers.get(i).add(sk);
                        kmers.get(i + 1).addAll(outKmers);

                        contigKmers.add(sk);
                    }

                    out.println("contigName: " + contigName);
                    out.println("sample: " + sample);
                    out.println("kmerSize: " + cg.getKmerSize());
                    out.println("orseq: " + StringUtil.repeatCharNTimes(' ', firstFlank.length()) + seq);
                    out.println("ieb02: " + new String(oeb[2]));
                    out.println("ieb01: " + new String(oeb[1]));
                    out.println("ieb00: " + new String(oeb[0]));
                    out.println("oedge: " + new String(oedges));
                    out.println("coseq: " + contig);
                    out.println("iedge: " + new String(iedges));
                    out.println("ieb00: " + new String(ieb[0]));
                    out.println("ieb01: " + new String(ieb[1]));
                    out.println("ieb02: " + new String(ieb[2]));
                    out.println("xaxis: " + new String(xaxis));

                    for (Integer pos : kmers.keySet()) {
                        for (String kmer : kmers.get(pos)) {
                            boolean inContig = contigKmers.contains(kmer);

                            CortexRecord cr = cg.findRecord(new CortexKmer(kmer));
                            int cov = (cr == null) ? 0 : cr.getCoverage(0);

                            out.println("kmer: " + Joiner.on(" ").join(pos, kmer, inContig, cov, linkMap.containsKey(kmer)));
                        }
                    }

                    List<String> iogt2l = new ArrayList<String>(iogt2);
                    for (int i = 0; i < iogt2l.size(); i++) {
                        for (int j = i + 1; j < iogt2l.size(); j++) {
                            String a = iogt2l.get(i);
                            String b = iogt2l.get(j);

                            if (!a.equals(b)) {
                                int cov = 0;

                                for (String sk : linkMap.keySet()) {
                                    CortexLinksRecord clr = linkMap.get(sk);
                                    for (CortexJunctionsRecord cjr : clr.getJunctions()) {
                                        Set<String> lks = kmersInLinks.get(cjr);

                                        if (lks.contains(a) && lks.contains(b)) {
                                            cov += cjr.getCoverage(0);
                                        }
                                    }
                                }

                                if (cov > 0) {
                                    out.println("pair: " + Joiner.on(" ").join(a, b, cov));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
