package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.assembly.cortex.CortexGraphWalker;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksMap;
import uk.ac.ox.well.indiana.utils.io.fasta.MultiIndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class PrintKnownVariantEdges extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference genome")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="bed", shortName="b", doc="Bed file")
    public File BED;

    /*
    @Argument(fullName="ref0", shortName="r0", doc="Reference genome for parent 0")
    public MultiIndexedFastaSequenceFile REF0;

    @Argument(fullName="ref1", shortName="r1", doc="Reference genome for parent 1")
    public MultiIndexedFastaSequenceFile REF1;
    */

    @Argument(fullName="graph", shortName="g", doc="Graph (sorted)")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Graph links")
    public File LINKS;

    @Argument(fullName="novelKmers", shortName="n", doc="Novel kmers")
    public FastaSequenceFile NOVEL_KMERS;

    @Output
    public PrintStream out;

    private class Link extends DefaultEdge {
        int coverage = 0;
    }

    private Map<String, String> parseInfoField(String info) {
        Map<String, String> infoMap = new HashMap<String, String>();
        String[] keyValuePairs = info.split(";");

        for (String keyValuePair : keyValuePairs) {
            String[] kv = keyValuePair.split("=");
            infoMap.put(kv[0], kv[1]);
        }

        return infoMap;
    }

    private Set<CortexKmer> loadNovelKmers() {
        Set<CortexKmer> novelKmers = new HashSet<CortexKmer>();
        ReferenceSequence rseq;
        while ((rseq = NOVEL_KMERS.nextSequence()) != null) {
            CortexKmer ck = new CortexKmer(rseq.getBases());

            novelKmers.add(ck);
        }

        return novelKmers;
    }

    private void printGraph(DirectedGraph<String, Link> g, PrintStream o) {
        o.println("digraph G {");

        for (String v : g.vertexSet()) {
            o.println("    \"" + v + "\";");
        }

        for (Link e : g.edgeSet()) {
            String s = g.getEdgeSource(e);
            String t = g.getEdgeTarget(e);

            o.println("    \"" + s + "\" -> \"" + t + "\" [ penwidth=" + Math.log10(e.coverage + 1) + " ];");
        }

        o.println("}");
    }

    private Set<String> walk(CortexGraph cg, CortexLinksMap linksMap, int color, String from, String to) {
        Set<String> reconstructionKmers = new HashSet<String>();
        StringBuilder reconstructionFwd = new StringBuilder();
        StringBuilder reconstructionRev = new StringBuilder();

        reconstructionKmers.add(from);
        reconstructionKmers.add(to);

        String thisKmer = from;
        reconstructionFwd.append(thisKmer);
        while (thisKmer != null) {
            Collection<String> kmers = CortexUtils.getPrevKmers(cg, thisKmer, color);

            if (kmers.size() == 1) {
                String tKmer = kmers.iterator().next();

                reconstructionKmers.add(tKmer);
                reconstructionFwd.insert(0, tKmer.charAt(0));

                thisKmer = tKmer;
            } else {
                for (String tKmer : kmers) { reconstructionKmers.add(tKmer); }
                thisKmer = null;
            }
        }

        thisKmer = from;
        while (thisKmer != null) {
            Collection<String> kmers = CortexUtils.getNextKmers(cg, thisKmer, color);

            if (kmers.size() == 1) {
                String tKmer = kmers.iterator().next();

                reconstructionKmers.add(tKmer);
                reconstructionFwd.append(tKmer.charAt(tKmer.length() - 1));

                thisKmer = tKmer;
            } else {
                for (String tKmer : kmers) { reconstructionKmers.add(tKmer); }
                thisKmer = null;
            }
        }

        thisKmer = to;
        reconstructionRev.append(thisKmer);
        while (thisKmer != null) {
            Collection<String> kmers = CortexUtils.getPrevKmers(cg, thisKmer, color);

            if (kmers.size() == 1) {
                String tKmer = kmers.iterator().next();

                reconstructionKmers.add(tKmer);
                reconstructionRev.insert(0, tKmer.charAt(0));

                thisKmer = tKmer;
            } else {
                for (String tKmer : kmers) { reconstructionKmers.add(tKmer); }
                thisKmer = null;
            }
        }

        thisKmer = to;
        while (thisKmer != null) {
            Collection<String> kmers = CortexUtils.getNextKmers(cg, thisKmer, color);

            if (kmers.size() == 1) {
                String tKmer = kmers.iterator().next();

                reconstructionKmers.add(tKmer);
                reconstructionRev.append(tKmer.charAt(tKmer.length() - 1));

                thisKmer = tKmer;
            } else {
                for (String tKmer : kmers) { reconstructionKmers.add(tKmer); }
                thisKmer = null;
            }
        }

        Set<String> alleles = new HashSet<String>();
        if (linksMap.getCortexLinks().hasColor(color)) {
            for (String reconstructionKmerFw : reconstructionKmers) {
                CortexKmer kmer = new CortexKmer(reconstructionKmerFw);
                if (linksMap.containsKey(kmer)) {
                    for (CortexJunctionsRecord cjr : linksMap.get(kmer).getJunctions()) {
                        String reconstruction = cjr.getSeq();

                        if (kmer.isFlipped() != cjr.isFlipped()) {
                            reconstruction = SequenceUtils.reverseComplement(reconstruction);
                        }

                        if (reconstruction.contains(from) && reconstruction.contains(to)) {
                            int fromIndex = reconstruction.indexOf(from) + from.length();
                            int toIndex = reconstruction.lastIndexOf(to);

                            String allele = reconstruction.substring(fromIndex, toIndex);
                            alleles.add(allele);
                        }
                    }
                }
            }
        } else {
            Set<String> reconstructions = new HashSet<String>();
            reconstructions.add(reconstructionFwd.toString());
            reconstructions.add(reconstructionRev.toString());

            for (String reconstruction : reconstructions) {
                if (reconstruction.contains(from) && reconstruction.contains(to)) {
                    int fromIndex = reconstruction.indexOf(from) + from.length();
                    int toIndex = reconstruction.lastIndexOf(to);

                    String allele = reconstruction.substring(fromIndex, toIndex);
                    alleles.add(allele);
                }
            }
        }

        return alleles;
    }

    private Set<String> walkRef(MultiIndexedFastaSequenceFile ref, String from, String to) {
        Interval fromInterval = ref.find(from);
        Interval toInterval = ref.find(to);

        Interval leftInterval  = fromInterval.getStart() < toInterval.getStart() ? fromInterval : toInterval;
        Interval rightInterval = fromInterval.getStart() < toInterval.getStart() ? toInterval   : fromInterval;

        Set<String> alleles = new HashSet<String>();
        if (fromInterval.getSequence().equals(toInterval.getSequence())) {
            String seq = new String(ref.getSubsequenceAt(leftInterval.getSequence(), leftInterval.getStart(), rightInterval.getEnd()).getBases());

            alleles.add(seq);
        }

        return alleles;
    }

    @Override
    public void execute() {
        Set<CortexKmer> novelKmers = loadNovelKmers();
        int kmerSize = novelKmers.iterator().next().length();

        CortexLinksMap linksMap = new CortexLinksMap(LINKS);

        TableReader tr = new TableReader(BED, "chr", "start", "stop", "info");

        Map<String, Integer> eventTypes = new HashMap<String, Integer>();
        int count = 0, full = 0;
        for (Map<String, String> te : tr) {
            String chr = te.get("chr");
            int start = Integer.valueOf(te.get("start")) + 1;
            int stop = Integer.valueOf(te.get("stop"));
            Map<String, String> infoMap = parseInfoField(te.get("info"));

            String eventType = infoMap.get("denovo");
            ContainerUtils.increment(eventTypes, eventType);

            if (!eventType.equals("unknown")) {
                String seq = new String(REF.getSubsequenceAt(chr, start, stop).getBases());
                String alt = infoMap.get("alt");

                int window = 200;

                int contextStart = (start - window >= 1) ? start - window : 1;
                int contextEnd = stop + window < REF.getSequence(chr).length() ? stop + window : REF.getSequence(chr).length() - 1;
                String context = new String(REF.getSubsequenceAt(chr, contextStart, contextEnd).getBases());
                StringBuilder novelty = new StringBuilder();
                StringBuilder edgesOut = new StringBuilder(StringUtil.repeatCharNTimes(' ', contextEnd - contextStart + kmerSize));
                StringBuilder edgesIn = new StringBuilder(StringUtil.repeatCharNTimes(' ', contextEnd - contextStart + kmerSize));

                for (int i = 0; i <= context.length() - kmerSize; i++) {
                    CortexKmer kmer = new CortexKmer(context.substring(i, i + kmerSize));
                    if (novelKmers.contains(kmer)) {
                        novelty.append(".");
                    } else {
                        novelty.append(" ");
                    }

                    CortexRecord cr = GRAPH.findRecord(kmer);
                    if (!kmer.isFlipped()) {
                        if (cr.getOutEdgesAsBytes(0).size() > 1) { edgesOut.setCharAt(i + kmerSize, '/'); }
                        else { edgesOut.setCharAt(i + kmerSize, ' '); }

                        if (cr.getInEdgesAsBytes(0).size() > 1) { edgesIn.setCharAt(i, '/'); }
                        else { edgesIn.setCharAt(i, ' '); }
                    } else {
                        if (cr.getOutEdgesAsBytes(0).size() > 1) { edgesIn.setCharAt(i, '/'); }
                        else { edgesIn.setCharAt(i, ' '); }

                        if (cr.getInEdgesAsBytes(0).size() > 1) { edgesOut.setCharAt(i + kmerSize, '/'); }
                        else { edgesOut.setCharAt(i + kmerSize, ' '); }
                    }
                }

                Interval leftInterval0 = null, leftInterval1 = null, rightInterval0 = null, rightInterval1 = null;

                for (int i = 0; i < window - kmerSize; i++) {
                    String kmer = context.substring(i, i + kmerSize);

                    if (leftInterval0 == null) {
                        Interval i0 = REF0.find(kmer);
                        if (i0 != null) {
                            leftInterval0 = i0;
                        }
                    }

                    if (leftInterval1 == null) {
                        Interval i1 = REF1.find(kmer);
                        if (i1 != null) {
                            leftInterval1 = i1;
                        }
                    }
                }

                for (int i = context.length() - kmerSize; i > context.length() - window; i--) {
                    String kmer = context.substring(i, i + kmerSize);

                    if (leftInterval0 != null && rightInterval0 == null) {
                        Interval i0 = REF0.find(kmer);
                        if (i0 != null && i0.getSequence().equals(leftInterval0.getSequence())) {
                            rightInterval0 = i0;
                        }
                    }

                    if (leftInterval1 != null && rightInterval1 == null) {
                        Interval i1 = REF1.find(kmer);
                        if (i1 != null && i1.getSequence().equals(leftInterval1.getSequence())) {
                            rightInterval1 = i1;
                        }
                    }
                }

                //String seq0 = (leftInterval0 == null || rightInterval0 == null) ? "" : new String(REF0.getSubsequenceAt(leftInterval0.getSequence(), leftInterval0.getStart(), rightInterval0.getEnd()).getBases());
                //String seq1 = (leftInterval1 == null || rightInterval1 == null) ? "" : new String(REF1.getSubsequenceAt(leftInterval1.getSequence(), leftInterval1.getStart(), rightInterval1.getEnd()).getBases());

                int leftFlankStart = start - kmerSize;
                int leftFlankEnd = start - 1;
                String leftFlank = new String(REF.getSubsequenceAt(chr, leftFlankStart, leftFlankEnd).getBases());

                int rightFlankStart = stop + 1;
                int rightFlankEnd = stop + kmerSize;
                String rightFlank = new String(REF.getSubsequenceAt(chr, rightFlankStart, rightFlankEnd).getBases());

                List<String> reconstructionKmers = new ArrayList<String>();

                String seq0 = "";
                if (leftInterval0 != null && rightInterval0 != null) {
                    if (leftInterval0.getStart() < rightInterval0.getStart()) {
                        seq0 = new String(REF0.getSubsequenceAt(leftInterval0.getSequence(), leftInterval0.getStart(), rightInterval0.getEnd()).getBases());
                    } else {
                        seq0 = new String(REF0.getSubsequenceAt(leftInterval0.getSequence(), rightInterval0.getStart(), leftInterval0.getEnd()).getBases());
                    }
                }

                if (seq0.equals("")) {
                    //walk(GRAPH, linksMap, 1, l)
                }

                String seq1 = "";
                if (leftInterval1 != null && rightInterval1 != null) {
                    if (leftInterval1.getStart() < rightInterval1.getStart()) {
                        seq1 = new String(REF1.getSubsequenceAt(leftInterval1.getSequence(), leftInterval1.getStart(), rightInterval1.getEnd()).getBases());
                    } else {
                        seq1 = new String(REF1.getSubsequenceAt(leftInterval1.getSequence(), rightInterval1.getStart(), leftInterval1.getEnd()).getBases());
                    }
                }

                /*
                String thisKmer = leftFlank;
                while (thisKmer != null) {
                    String nextKmer = CortexUtils.getNextKmer(GRAPH, thisKmer);

                    if (nextKmer != null) {
                        reconstructionKmers.add(nextKmer);
                    }

                    thisKmer = nextKmer;
                }

                StringBuilder reconstruction = new StringBuilder();

                if (reconstructionKmers.size() >= 2) {
                    for (int i = 0; i < reconstructionKmers.size(); i++) {
                        reconstruction.append(reconstructionKmers.get(i).charAt(kmerSize - 1));
                    }
                }

                Set<String> trialFlanks = new HashSet<String>(Arrays.asList(leftFlank));
                String commonLeftFlankRef0 = leftFlank;
                while (trialFlanks.size() == 1 && GRAPH.findRecord(new CortexKmer(commonLeftFlankRef0)).getCoverage(1) != 1) {
                    trialFlanks = CortexUtils.getPrevKmers(GRAPH, commonLeftFlankRef0, 0);

                    if (trialFlanks.size() == 1) {
                        commonLeftFlankRef0 = trialFlanks.iterator().next();
                    }
                }

                trialFlanks = new HashSet<String>(Arrays.asList(rightFlank));
                String commonRightFlankRef0 = rightFlank;
                while (trialFlanks.size() == 1 && GRAPH.findRecord(new CortexKmer(commonRightFlankRef0)).getCoverage(1) != 1) {
                    trialFlanks = CortexUtils.getNextKmers(GRAPH, commonRightFlankRef0, 0);

                    if (trialFlanks.size() == 1) {
                        commonRightFlankRef0 = trialFlanks.iterator().next();
                    }
                }

                trialFlanks = new HashSet<String>(Arrays.asList(leftFlank));
                String commonLeftFlankRef1 = leftFlank;
                while (trialFlanks.size() == 1 && GRAPH.findRecord(new CortexKmer(commonLeftFlankRef1)).getCoverage(2) != 1) {
                    trialFlanks = CortexUtils.getPrevKmers(GRAPH, commonLeftFlankRef0, 0);

                    if (trialFlanks.size() == 1) {
                        commonLeftFlankRef1 = trialFlanks.iterator().next();
                    }
                }

                trialFlanks = new HashSet<String>(Arrays.asList(rightFlank));
                String commonRightFlankRef1 = rightFlank;
                while (trialFlanks.size() == 1 && GRAPH.findRecord(new CortexKmer(commonRightFlankRef1)).getCoverage(2) != 1) {
                    trialFlanks = CortexUtils.getNextKmers(GRAPH, commonRightFlankRef1, 0);

                    if (trialFlanks.size() == 1) {
                        commonRightFlankRef1 = trialFlanks.iterator().next();
                    }
                }
                */

                log.debug("event: {}",   eventType);
                log.debug("  seq: {}{}", StringUtil.repeatCharNTimes(' ', start - contextStart), seq);
                log.debug("  alt: {}{}", StringUtil.repeatCharNTimes(' ', start - contextStart), alt);
                log.debug("  out: {}{}", StringUtil.repeatCharNTimes(' ', kmerSize), edgesOut);
                log.debug("  con: {}",   context);
                log.debug("   in: {}",   edgesIn);
                log.debug("  nov: {}",   novelty);
                log.debug("   s0: {}",   seq0);
                log.debug("   s1: {}",   seq1);

//                log.info("   lf: {}{}", StringUtil.repeatCharNTimes(' ', leftFlankStart - contextStart), leftFlank);
//                log.info("   rf: {}{}", StringUtil.repeatCharNTimes(' ', rightFlankStart - contextStart), rightFlank);
//                log.info("  lf0: {} {}", commonLeftFlankRef0, GRAPH.findRecord(new CortexKmer(commonLeftFlankRef0)));
//                log.info("  lf1: {} {}", commonLeftFlankRef1, GRAPH.findRecord(new CortexKmer(commonLeftFlankRef1)));
//                log.info("  rf0: {} {}", commonRightFlankRef0, GRAPH.findRecord(new CortexKmer(commonRightFlankRef0)));
//                log.info("  rf1: {} {}", commonRightFlankRef1, GRAPH.findRecord(new CortexKmer(commonRightFlankRef1)));

                /*
                DirectedGraph<String, Link> g = new DefaultDirectedGraph<String, Link>(Link.class);

                for (int i = 0; i <= context.length() - kmerSize; i++) {
                    String tk = context.substring(i, i + kmerSize);
                    g.addVertex(tk);
                    if (i < context.length() - kmerSize) {
                        String nk = context.substring(i + 1, i + 1 + kmerSize);
                        g.addVertex(nk);
                        g.addEdge(tk, nk, new Link());
                    }

                    CortexKmer kmer = new CortexKmer(context.substring(i, i + kmerSize));
                    if (linksMap.containsKey(kmer)) {
                        for (CortexJunctionsRecord cjr : linksMap.get(kmer).getJunctions()) {
                            String linkSeq = cjr.getSeq();
                            int padding = i;

                            if (kmer.isFlipped() != cjr.isFlipped()) {
                                linkSeq = SequenceUtils.reverseComplement(linkSeq);
                                padding = i - linkSeq.length() + kmerSize;
                            }

                            StringBuilder sb = new StringBuilder(linkSeq);
                            if (padding < 0) {
                                sb.delete(0, -padding);
                                padding = 0;
                            }

                            linkSeq = sb.toString();

                            if (linkSeq.contains(leftFlank) || linkSeq.contains(rightFlank)) {
                                log.info("  lnk: {}{}", StringUtil.repeatCharNTimes(' ', padding), sb);
                            }
                        }
                    }
                }
                */

                /*
                Set<String> childAlleles0 = walk(GRAPH, linksMap, 0, commonLeftFlankRef0, commonRightFlankRef0);
                Set<String> ref0Alleles = walkRef(REF0, commonLeftFlankRef0, commonRightFlankRef0);
                Set<String> childAlleles1 = walk(GRAPH, linksMap, 0, commonLeftFlankRef1, commonRightFlankRef1);
                Set<String> ref1Alleles = walkRef(REF1, commonLeftFlankRef1, commonRightFlankRef1);

                log.info("child: {}", Joiner.on(", ").join(childAlleles0));
                log.info(" ref0: {}", Joiner.on(", ").join(ref0Alleles));
                log.info("child: {}", Joiner.on(", ").join(childAlleles1));
                log.info(" ref1: {}", Joiner.on(", ").join(ref1Alleles));

                log.info("--");
                */
                count++;
                if (!seq0.equals("") && !seq1.equals("")) {
                    full++;
                }

                log.info("{} {} {}:{}-{} {}", count, full, chr, start, stop, eventType);
            }
        }

        log.info("  {}", Joiner.on("; ").withKeyValueSeparator("=").join(eventTypes));
    }
}
