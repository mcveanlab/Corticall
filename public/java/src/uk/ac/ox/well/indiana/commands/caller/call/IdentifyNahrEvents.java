package uk.ac.ox.well.indiana.commands.caller.call;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;

import java.io.PrintStream;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by kiran on 22/06/2017.
 */
public class IdentifyNahrEvents extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="refs", shortName="R", doc="References")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Argument(fullName="sequences", shortName="s", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="validated", shortName="v", doc="Validated NAHR event")
    public FastaSequenceFile NAHR;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(PARENTS);
        List<Integer> recruitColors = GRAPH.getColorsForSampleNames(new ArrayList<>(LOOKUPS.keySet()));

        Set<Integer> colors = new TreeSet<>();
        colors.add(childColor);
        colors.addAll(parentColors);
        colors.addAll(recruitColors);

        Map<CortexKmer, Set<String>> expectedNovelKmers = new HashMap<>();
        ReferenceSequence rseq;
        while ((rseq = NAHR.nextSequence()) != null) {
            String seq = rseq.getBaseString();
            for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                CortexKmer ck = new CortexKmer(seq.substring(i, i + GRAPH.getKmerSize()));
                CortexRecord rr = ROI.findRecord(ck);

                if (rr != null) {
                    expectedNovelKmers.put(rr.getCortexKmer(), new TreeSet<>());
                }
            }
        }

        Map<String, String> contigs = new HashMap<>();

        while ((rseq = CONTIGS.nextSequence()) != null) {
            String[] name = rseq.getName().split("\\s+");
            String seq = rseq.getBaseString();
            contigs.put(name[0], seq);

            for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                CortexKmer ck = new CortexKmer(seq.substring(i, i + GRAPH.getKmerSize()));

                if (expectedNovelKmers.containsKey(ck)) {
                    ContainerUtils.add(expectedNovelKmers, ck, name[0]);
                }
            }
        }

        Set<String> cnames = new HashSet<>();
        for (CortexKmer ck : expectedNovelKmers.keySet()) {
            cnames.addAll(expectedNovelKmers.get(ck));
        }

        Map<String, String> enc = createContigEncoding();

        Set<ContigInfo> candidates = findCandidateNAHRs(contigs, colors, enc);
        Map<CortexKmer, ContigInfo> cEnds = new HashMap<>();

        for (ContigInfo ci : candidates) {
            cEnds.put(ci.crFirst.getCortexKmer(), ci);
            cEnds.put(ci.crLast.getCortexKmer(), ci);
        }

        Set<String> usedCandidates = new HashSet<>();

        String firstChrPattern = "([^_\\.])";
        Pattern firstChrMotif = Pattern.compile(firstChrPattern);

        String lastChrPattern = "([^_\\.])[_\\.]*$";
        Pattern lastChrMotif = Pattern.compile(lastChrPattern);

        for (ContigInfo candidate : candidates) {
            if (!usedCandidates.contains(candidate.name)) {
                usedCandidates.add(candidate.name);

                Matcher firstChrMatcher = firstChrMotif.matcher(candidate.anntig);

                List<String> bridgeBack = new ArrayList<>();
                String bridgeBackName = "";
                if (firstChrMatcher.find()) {
                    String encchr = firstChrMatcher.group(1);
                    String chrName = null;
                    for (String n : enc.keySet()) {
                        if (enc.get(n).equals(encchr)) {
                            chrName = n;
                        }
                    }

                    if (chrName != null) {
                        String sk = candidate.skFirst;
                        do {
                            CortexRecord cr = GRAPH.findRecord(new CortexKmer(sk));
                            Map<Integer, Set<String>> pks = TraversalEngine.getAllPrevKmers(cr, !sk.equals(cr.getKmerAsString()));
                            String pk = null;
                            for (String p : pks.get(childColor)) {
                                for (String background : LOOKUPS.keySet()) {
                                    if (LOOKUPS.get(background).getReferenceSequence().getSequenceDictionary().getSequence(chrName) != null) {
                                        Set<Interval> its = LOOKUPS.get(background).findKmer(p);

                                        if (its.size() == 1) {
                                            Interval it = its.iterator().next();
                                            if (it.getContig().equals(chrName)) {
                                                pk = p;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }

                            sk = null;

                            if (pk != null) {
                                CortexKmer ck = new CortexKmer(pk);
                                if (cEnds.containsKey(ck)) {
                                    ContigInfo cci = cEnds.get(ck);

                                    if (pk.equals(cci.skLast)) {
                                        bridgeBack.add(0, cci.contig);
                                    } else if (pk.equals(SequenceUtils.reverseComplement(cci.skFirst))) {
                                        bridgeBack.add(0, SequenceUtils.reverseComplement(cci.contig));
                                    }

                                    usedCandidates.add(cci.name);
                                    bridgeBackName = cci.name;
                                } else {
                                    bridgeBack.add(0, pk);

                                    sk = pk;
                                }
                            }
                        } while (sk != null);
                    }
                }

                Matcher lastChrMatcher = lastChrMotif.matcher(candidate.anntig);
                List<String> bridgeFwd = new ArrayList<>();
                String bridgeFwdName = "";
                if (lastChrMatcher.find()) {
                    String encchr = lastChrMatcher.group(1);
                    String chrName = null;
                    for (String n : enc.keySet()) {
                        if (enc.get(n).equals(encchr)) {
                            chrName = n;
                        }
                    }

                    if (chrName != null) {
                        String sk = candidate.skLast;
                        do {
                            CortexRecord cr = GRAPH.findRecord(new CortexKmer(sk));
                            Map<Integer, Set<String>> nks = TraversalEngine.getAllNextKmers(cr, !sk.equals(cr.getKmerAsString()));
                            String nk = null;
                            for (String n : nks.get(childColor)) {
                                for (String background : LOOKUPS.keySet()) {
                                    if (LOOKUPS.get(background).getReferenceSequence().getSequenceDictionary().getSequence(chrName) != null) {
                                        Set<Interval> its = LOOKUPS.get(background).findKmer(n);

                                        if (its.size() == 1) {
                                            Interval it = its.iterator().next();
                                            if (it.getContig().equals(chrName)) {
                                                nk = n;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }

                            sk = null;

                            if (nk != null) {
                                CortexKmer ck = new CortexKmer(nk);
                                if (cEnds.containsKey(ck)) {
                                    ContigInfo cci = cEnds.get(ck);

                                    if (nk.equals(cci.skFirst)) {
                                        bridgeFwd.add(cci.contig);
                                    } else if (nk.equals(SequenceUtils.reverseComplement(cci.skLast))) {
                                        bridgeFwd.add(SequenceUtils.reverseComplement(cci.contig));
                                    }

                                    usedCandidates.add(cci.name);
                                    bridgeFwdName = cci.name;
                                } else {
                                    bridgeBack.add(nk);

                                    sk = nk;
                                }
                            }
                        } while (sk != null);
                    }
                }

                if (bridgeBack.size() > 0 || bridgeFwd.size() > 0) {
                    StringBuilder newContig = new StringBuilder(candidate.contig);
                    StringBuilder prev = new StringBuilder();
                    for (String p : bridgeBack) {
                        if (prev.length() == 0) {
                            prev.append(p);
                        } else {
                            prev.append(p.substring(p.length() - 1, p.length()));
                        }
                    }

                    StringBuilder next = new StringBuilder();
                    for (int i = 0; i < bridgeFwd.size(); i++) {
                        String n = bridgeFwd.get(i);
                        if (i < bridgeFwd.size() - 1) {
                            next.append(n.substring(0, 1));
                        } else {
                            next.append(n);
                        }
                    }

                    newContig.insert(0, prev.toString());
                    newContig.append(next.toString());

                    log.info("joined: {} {} {} {}", candidate.contig.length(), newContig.length(), bridgeBackName, bridgeFwdName);
                }
            }
        }
    }

    private class ContigInfo {
        public String name;
        public String contig;
        public String anntig;
        public String skFirst;
        public String skLast;
        public CortexRecord crFirst;
        public CortexRecord crLast;

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            ContigInfo that = (ContigInfo) o;

            if (name != null ? !name.equals(that.name) : that.name != null) return false;
            if (contig != null ? !contig.equals(that.contig) : that.contig != null) return false;
            if (anntig != null ? !anntig.equals(that.anntig) : that.anntig != null) return false;
            if (crFirst != null ? !crFirst.equals(that.crFirst) : that.crFirst != null) return false;
            return crLast != null ? crLast.equals(that.crLast) : that.crLast == null;

        }

        @Override
        public int hashCode() {
            int result = name != null ? name.hashCode() : 0;
            result = 31 * result + (contig != null ? contig.hashCode() : 0);
            result = 31 * result + (anntig != null ? anntig.hashCode() : 0);
            result = 31 * result + (crFirst != null ? crFirst.hashCode() : 0);
            result = 31 * result + (crLast != null ? crLast.hashCode() : 0);
            return result;
        }
    }

    private Set<ContigInfo> findCandidateNAHRs(Map<String, String> contigs, Set<Integer> colors, Map<String, String> enc) {
        String flankingNovelPattern = "(([^_\\.])\\2+)_*(\\.+)_*(([^_\\.])\\5+)";
        Pattern flankingNovelMotif = Pattern.compile(flankingNovelPattern);

        String novelPattern = "(\\.+)";
        Pattern novelMotif = Pattern.compile(novelPattern);

        Set<ContigInfo> candidates = new HashSet<>();
        for (String cname : contigs.keySet()) {
            String contig = contigs.get(cname);

            for (String background : LOOKUPS.keySet()) {
                String anntig = annotateContig(LOOKUPS.get(background), contig, enc);

                Matcher flankingNovelMatcher = flankingNovelMotif.matcher(anntig);
                int numRecombs = 0;
                if (flankingNovelMatcher.find()) {
                    do {
                        if (!flankingNovelMatcher.group(2).equals(flankingNovelMatcher.group(5))) {
                            numRecombs++;
                        }
                    } while (flankingNovelMatcher.find(flankingNovelMatcher.start(3)));
                }

                Matcher novelMatcher = novelMotif.matcher(anntig);
                int numNovelRuns = 0;
                while (novelMatcher.find()) {
                    numNovelRuns++;
                }

                if (numRecombs > 0 && numNovelRuns > 1) {
                    String skFirst = contig.substring(0, GRAPH.getKmerSize());
                    CortexKmer ckFirst = new CortexKmer(skFirst);
                    CortexRecord crFirst = GRAPH.findRecord(ckFirst);

                    String skLast = contig.substring(contig.length() - GRAPH.getKmerSize(), contig.length());
                    CortexKmer ckLast = new CortexKmer(skLast);
                    CortexRecord crLast = GRAPH.findRecord(ckLast);

                    log.info("name {} length {}", cname, contig.length());
                    log.info("    numRecombs={} numNovelRuns={}", numRecombs, numNovelRuns);
                    log.info("    contig     {}", contig);
                    log.info("    anntig {} {}", background, anntig);
                    log.info("    skFirst {}", skFirst);
                    log.info("    skLast  {}", skLast);
                    log.info("    crFirst {}", recordToString(skFirst, crFirst, colors));
                    log.info("    crLast  {}", recordToString(skLast, crLast, colors));

                    ContigInfo ci = new ContigInfo();
                    ci.name = cname;
                    ci.contig = contig;
                    ci.anntig = anntig;
                    ci.skFirst = skFirst;
                    ci.skLast = skLast;
                    ci.crFirst = crFirst;
                    ci.crLast = crLast;
                    candidates.add(ci);
                }
            }
        }

        return candidates;
    }

    private String annotateContig(KmerLookup kl, String contig, Map<String, String> contigEncoding) {
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i <= contig.length() - GRAPH.getKmerSize(); i++) {
            String sk = contig.substring(i, i + GRAPH.getKmerSize());
            CortexKmer ck = new CortexKmer(sk);
            CortexRecord rr = ROI.findRecord(ck);
            Set<Interval> loci = kl.findKmer(sk);

            if (rr != null) {
                sb.append(".");
            } else if (loci.size() == 1) {
                Interval locus = loci.iterator().next();
                sb.append(contigEncoding.get(locus.getContig()));
            } else {
                sb.append("_");
            }
        }

        return sb.toString();
    }

    @NotNull
    private Map<String, String> createContigEncoding() {
        Map<String, String> contigEncoding = new HashMap<>();
        String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
        Random r = new Random(0);
        for (String key : LOOKUPS.keySet()) {
            Set<String> usedCodes = new HashSet<>();

            for (SAMSequenceRecord ssr : LOOKUPS.get(key).getReferenceSequence().getSequenceDictionary().getSequences()) {
                String sn = ssr.getSequenceName();
                String c;
                do {
                    c = String.valueOf(alphabet.charAt(r.nextInt(alphabet.length())));
                } while(usedCodes.contains(c) && usedCodes.size() < alphabet.length());

                if (!usedCodes.contains(c)) {
                    contigEncoding.put(sn, c);
                    usedCodes.add(c);
                }
            }
        }
        return contigEncoding;
    }

    private String recordToString(String sk, CortexRecord cr, Set<Integer> colors) {
        String kmer = cr.getKmerAsString();
        String cov = "";
        String ed = "";

        boolean fw = sk.equals(kmer);

        if (!fw) {
            kmer = SequenceUtils.reverseComplement(kmer);
        }

        int color = 0;
        for (int coverage : cr.getCoverages()) {
            if (colors.contains(color)) {
                cov += " " + coverage;
            }
            color++;
        }

        color = 0;
        for (String edge : cr.getEdgeAsStrings()) {
            if (colors.contains(color)) {
                ed += " " + (fw ? edge : SequenceUtils.reverseComplement(edge));
            }
            color++;
        }

        /*
        Set<String> lss = new TreeSet<>();
        if (LOOKUPS != null) {
            Set<Interval> loci = LOOKUP.findKmer(kmer);

            if (loci != null && loci.size() > 0) {
                for (Interval locus : loci) {
                    String ls = locus.getContig() + ":" + locus.getStart() + "-" + locus.getEnd() + ":" + (locus.isPositiveStrand() ? "+" : "-");
                    lss.add(ls);
                }
            }
        }
        String lssCombined = Joiner.on(";").join(lss);

        return kmer + " " + cov + " " + ed + " " + lssCombined;
        */

        return kmer + " " + cov + " " + ed;
    }
}
