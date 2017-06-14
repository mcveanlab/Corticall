package uk.ac.ox.well.indiana.commands.caller.call;

import com.google.api.client.util.Joiner;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.jgrapht.traverse.TopologicalOrderIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.stoppingconditions.NahrStopper;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;
import uk.ac.ox.well.indiana.utils.visualizer.GraphVisualizer;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

/**
 * Created by kiran on 05/06/2017.
 */
public class CallNAHRs extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="lookups", shortName="l", doc="Lookups")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Output(fullName="out0", shortName="o0", doc="Out 0")
    public File f0;

    @Output(fullName="out1", shortName="o1", doc="Out 1")
    public File f1;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(PARENTS);

        SAMFileHeader sfh0 = new SAMFileHeader();
        sfh0.setSequenceDictionary(LOOKUPS.get("ref").getReferenceSequence().getSequenceDictionary());
        sfh0.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        SAMFileWriter sfw0 = new SAMFileWriterFactory().makeBAMWriter(sfh0, false, f0);

        SAMFileHeader sfh1 = new SAMFileHeader();
        sfh1.setSequenceDictionary(LOOKUPS.get("HB3").getReferenceSequence().getSequenceDictionary());
        sfh1.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        SAMFileWriter sfw1 = new SAMFileWriterFactory().makeBAMWriter(sfh1, false, f1);

        Map<CortexKmer, Boolean> used = new HashMap<>();
        for (CortexRecord rr : ROI) {
            used.put(rr.getCortexKmer(), false);
        }

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing novel kmers")
                .message("kmers processed")
                .maxRecord(used.size())
                .updateRecord(used.size())
                .make(log);

        for (CortexKmer ck : used.keySet()) {
            if (!used.get(ck)) {
                String sk = ck.getKmerAsString();

                Pair<List<String>, List<Interval>> recon0 = reconstruct("ref", sk);
                Pair<List<String>, List<Interval>> recon1 = reconstruct("HB3", sk);

                Set<Interval> mergedIntervals0 = mergeIntervals(recon0);
                Set<Interval> mergedIntervals1 = mergeIntervals(recon1);

                Map<String, Interval> aggregatedIntervals0 = aggregateIntervals(mergeIntervals(recon0));
                Map<String, Interval> aggregatedIntervals1 = aggregateIntervals(mergeIntervals(recon1));

                int novels0 = 0;
                for (String qk : recon0.getFirst()) {
                    CortexKmer rk = new CortexKmer(qk);
                    if (used.containsKey(rk)) { used.put(rk, true); novels0++; }
                }

                int novels1 = 0;
                for (String qk : recon1.getFirst()) {
                    CortexKmer rk = new CortexKmer(qk);
                    if (used.containsKey(rk)) { used.put(rk, true); novels1++; }
                }

                for (int i = 0; i < Math.max(novels0, novels1); i++) {
                    pm.update();
                }

                if (novels0 >= 10 && aggregatedIntervals0.size() > 1 && hasMultiChrBreakpoint(recon0, used)) {
                    List<SAMRecord> contig0 = getContig(recon0, sfh0, sk);

                    for (SAMRecord r : contig0) { sfw0.addAlignment(r); }

                    log.info("  {}", sk);
                    log.info("    - {} {} {} {}", recon0.getFirst().size(), novels0, mergedIntervals0);
                    log.info("      {}", contig0);

                    //printReconstruction(recon0, aggregatedIntervals0, "ref");
                }

                if (novels1 >= 10 && aggregatedIntervals1.size() > 1 && hasMultiChrBreakpoint(recon1, used)) {
                    List<SAMRecord> contig1 = getContig(recon1, sfh1, sk);

                    for (SAMRecord r : contig1) { sfw1.addAlignment(r); }

                    log.info("  {}", sk);
                    log.info("    - {} {} {} {}", recon1.getFirst().size(), novels1, mergedIntervals1);
                    log.info("      {}", contig1);

                    //printReconstruction(recon1, aggregatedIntervals1, "HB3");
                }
            }
        }

        sfw0.close();
        sfw1.close();
    }

    private List<SAMRecord> getContig(Pair<List<String>, List<Interval>> recon, SAMFileHeader sfh, String kmer) {
        List<SAMRecord> srs = new ArrayList<>();

        StringBuilder sb = new StringBuilder();
        Interval it0 = null;
        int numNovels = 0;

        for (int i = 0; i < recon.getFirst().size(); i++) {
            String sk = recon.getFirst().get(i);
            Interval it1 = recon.getSecond().get(i);

            sb.append(sb.length() == 0 ? sk : sk.substring(sk.length() - 1, sk.length()));

            if (it1 != null) {
                if (it0 == null) {
                    it0 = it1;
                } else if (it0.getIntersectionLength(it1) == GRAPH.getKmerSize() - 1) {
                    it0 = new Interval(it0.getContig(), it0.getStart() < it1.getStart() ? it0.getStart() : it1.getStart(), it0.getEnd() > it1.getEnd() ? it0.getEnd() : it1.getEnd(), it0.isNegativeStrand(), null);
                } else {
                    List<CigarElement> ces = new ArrayList<>();
                    ces.add(new CigarElement(sb.length() - numNovels, CigarOperator.M));
                    if (numNovels > 0) {
                        ces.add(new CigarElement(numNovels, CigarOperator.S));
                    }
                    Cigar cigar = new Cigar(ces);

                    SAMRecord sr = new SAMRecord(sfh);
                    sr.setReadName("read_" + kmer);
                    sr.setReadBases(sb.toString().getBytes());
                    sr.setReferenceName(it0.getContig());
                    sr.setAlignmentStart(it0.getStart());
                    sr.setCigar(cigar);
                    sr.setReadNegativeStrandFlag(it0.isNegativeStrand());
                    sr.setHeader(sfh);
                    srs.add(sr);

                    sb = new StringBuilder();
                    it0 = null;
                }
            } else {
                numNovels++;
            }
        }

        if (it0 != null) {
            List<CigarElement> ces = new ArrayList<>();
            ces.add(new CigarElement(sb.length() - numNovels, CigarOperator.M));
            if (numNovels > 0) {
                ces.add(new CigarElement(numNovels, CigarOperator.S));
            }
            Cigar cigar = new Cigar(ces);

            SAMRecord sr = new SAMRecord(sfh);
            sr.setReadName("read_" + kmer);
            sr.setReadBases(sb.toString().getBytes());
            sr.setReferenceName(it0.getContig());
            sr.setAlignmentStart(it0.getStart());
            sr.setCigar(cigar);
            sr.setReadNegativeStrandFlag(it0.isNegativeStrand());
            sr.setHeader(sfh);
            srs.add(sr);
        }

        return srs;
    }

    private boolean hasMultiChrBreakpoint(Pair<List<String>, List<Interval>> recon, Map<CortexKmer, Boolean> used) {
        for (int i = 1; i < recon.getFirst().size() - 1; i++) {
            CortexKmer ck = new CortexKmer(recon.getFirst().get(i));

            if (used.containsKey(ck)) {
                Interval b0 = null;
                for (int j = i - 1; j >= 0; j--) {
                    if (recon.getSecond().get(j) != null) {
                        b0 = recon.getSecond().get(j);
                        break;
                    }
                }

                Interval b1 = null;
                for (int j = i + 1; j < recon.getFirst().size(); j++) {
                    if (recon.getSecond().get(j) != null) {
                        b1 = recon.getSecond().get(j);
                        i = j;
                        break;
                    }
                }

                if (b0 != null && b1 != null && !b0.getContig().equals(b1.getContig())) {
                    return true;
                }
            }
        }

        return false;
    }

    private void printReconstruction(Pair<List<String>, List<Interval>> recon, Map<String, Interval> aggregatedIntervals, String background) {
        Map<String, Integer> contigIndices = new HashMap<>();
        int index = aggregatedIntervals.size();
        for (String contig : aggregatedIntervals.keySet()) {
            contigIndices.put(contig, index);
            index--;
        }

        out.println(Joiner.on('\t').join(Arrays.asList("kmer", "interval", "pos", "contigIndex")));

        for (int i = 0; i < recon.getFirst().size(); i++) {
            String kmer = recon.getFirst().get(i);
            Interval interval = recon.getSecond().get(i);
            int contigIndex = interval == null ? -1 : contigIndices.get(interval.getContig());

            log.info("{} {} {} {}", i, kmer, interval, LOOKUPS.get(background).findKmer(kmer));

            if (contigIndex >= 0) {
                String intervalString = interval.getContig() + ":" + interval.getStart() + "-" + interval.getEnd() + ":" + (interval.isPositiveStrand() ? "+" : "-");
                out.println(Joiner.on('\t').join(Arrays.asList(kmer, intervalString, i, contigIndex)));
            } else {
                out.println(Joiner.on('\t').join(Arrays.asList(kmer, "NA", i, contigIndex)));
            }
        }
    }

    private DirectedWeightedPseudograph<CortexVertex, CortexEdge> buildGraph(Pair<List<String>, List<Interval>> recon) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = new DirectedWeightedPseudograph<>(CortexEdge.class);

        for (int i = 0; i < recon.getFirst().size() - 1; i++) {
            String s0 = recon.getFirst().get(i);
            CortexRecord c0 = GRAPH.findRecord(new CortexKmer(s0));
            Interval i0 = recon.getSecond().get(i);
            CortexVertex v0 = new CortexVertex(s0, c0, i0);

            String s1 = recon.getFirst().get(i + 1);
            CortexRecord c1 = GRAPH.findRecord(new CortexKmer(s1));
            Interval i1 = recon.getSecond().get(i + 1);
            CortexVertex v1 = new CortexVertex(s1, c1, i1);

            g.addVertex(v0);
            g.addVertex(v1);

            Map<Integer, Set<String>> aks = TraversalEngine.getAllNextKmers(c0, !s0.equals(c0.getKmerAsString()));
            Set<String> pnk = aks.get(GRAPH.getColorForSampleName("ref"));
            Set<String> cnk = aks.get(GRAPH.getColorForSampleName(CHILD));

            if (pnk.contains(s1)) {
                g.addEdge(v0, v1, new CortexEdge(GRAPH.getColorForSampleName("ref"), 1.0));
            }

            if (cnk.contains(s1)) {
                g.addEdge(v0, v1, new CortexEdge(GRAPH.getColorForSampleName(CHILD), 1.0));
            }
        }

        return g;
    }

    @NotNull
    private Map<String, Interval> aggregateIntervals(Set<Interval> mergedIntervals) {
        Map<String, Interval> aggregatedIntervals = new TreeMap<>();
        for (Interval it : mergedIntervals) {
            if (!aggregatedIntervals.containsKey(it.getContig())) {
                aggregatedIntervals.put(it.getContig(), it);
            } else {
                Interval locus = aggregatedIntervals.get(it.getContig());

                int start = locus.getStart() < it.getStart() ? locus.getStart() : it.getStart();
                int end   = locus.getEnd()   > it.getEnd()   ? locus.getEnd()   : it.getEnd();

                Interval newlocus = new Interval(locus.getContig(), start, end, locus.isNegativeStrand(), null);

                aggregatedIntervals.put(locus.getContig(), newlocus);
            }
        }
        return aggregatedIntervals;
    }

    @NotNull
    private Set<Interval> mergeIntervals(Pair<List<String>, List<Interval>> recon) {
        Set<Interval> mergedIntervals = new LinkedHashSet<>();

        Interval locus = null;
        for (Interval it : recon.getSecond()) {
            if (locus == null) {
                locus = it;
            } else {
                if (it != null && locus.getContig().equals(it.getContig()) && locus.isPositiveStrand() == it.isPositiveStrand() && locus.getIntersectionLength(it) == GRAPH.getKmerSize() - 1) {
                    int start = locus.getStart() < it.getStart() ? locus.getStart() : it.getStart();
                    int end   = locus.getEnd()   > it.getEnd()   ? locus.getEnd()   : it.getEnd();

                    locus = new Interval(locus.getContig(), start, end, locus.isNegativeStrand(), null);
                } else {
                    mergedIntervals.add(locus);

                    locus = it;
                }
            }
        }

        if (locus != null) {
            mergedIntervals.add(locus);
        }

        return mergedIntervals;
    }

    private Pair<List<String>, List<Interval>> reconstruct(String background, String sk) {
        Pair<List<String>, List<Interval>> rev = reconstruct(background, sk, false, 4000);
        Pair<List<String>, List<Interval>> fwd = reconstruct(background, sk, true, 4000);

        List<String> allKmers = new ArrayList<>();
        List<Interval> allLoci = new ArrayList<>();

        allKmers.addAll(rev.getFirst());
        allKmers.add(sk);
        allKmers.addAll(fwd.getFirst());

        allLoci.addAll(rev.getSecond());
        allLoci.add(null);
        allLoci.addAll(fwd.getSecond());

        return new Pair<>(allKmers, allLoci);
    }

    private Pair<List<String>, List<Interval>> reconstruct(String background, String sk, boolean goForward, int limit) {
        List<String> vertices = new ArrayList<>();
        List<Interval> loci = new ArrayList<>();

        Set<String> usedNovelKmers = new HashSet<>();
        Set<Interval> usedLoci = new HashSet<>();

        Interval ci = null;
        int distanceFromNovel = 0;
        boolean onRef = false;
        boolean positiveStrand = false;
        boolean keepGoing = true;

        log.debug("start");

        while (distanceFromNovel < limit && keepGoing) {
            if (ci == null) {
                log.debug("  {} {} {} {}", goForward, distanceFromNovel, sk, ci);

                CortexRecord cr = GRAPH.findRecord(new CortexKmer(sk));
                Map<Integer, Set<String>> aks = goForward ? TraversalEngine.getAllNextKmers(cr, !sk.equals(cr.getKmerAsString())) : TraversalEngine.getAllPrevKmers(cr, !sk.equals(cr.getKmerAsString()));
                Set<String> achi = aks.get(GRAPH.getColorForSampleName(CHILD));

                if (achi.size() == 1) {
                    sk = achi.iterator().next();

                    Set<Interval> acis = LOOKUPS.get(background).findKmer(sk);
                    ci = acis.size() == 1 ? acis.iterator().next() : null;

                    if (goForward) {
                        vertices.add(sk);
                        loci.add(ci);
                    } else {
                        vertices.add(0, sk);
                        loci.add(0, ci);
                    }

                    if (ci != null) {
                        onRef = true;
                        positiveStrand = ci.isPositiveStrand();
                    } else if (!usedNovelKmers.contains(sk)) {
                        distanceFromNovel = 0;
                        usedNovelKmers.add(sk);
                    } else {
                        keepGoing = false;
                    }
                } else {
                    keepGoing = false;
                }
            } else {
                log.debug("  {} {} {} {}", goForward, distanceFromNovel, sk, ci);

                do {
                    Interval aci;

                    if (positiveStrand) {
                        if (goForward) {
                            aci = new Interval(ci.getContig(), ci.getStart() + 1, ci.getEnd() + 1, ci.isNegativeStrand(), null);
                        } else {
                            aci = new Interval(ci.getContig(), ci.getStart() - 1, ci.getEnd() - 1, ci.isNegativeStrand(), null);
                        }
                    } else {
                        if (goForward) {
                            aci = new Interval(ci.getContig(), ci.getStart() - 1, ci.getEnd() - 1, ci.isNegativeStrand(), null);
                        } else {
                            aci = new Interval(ci.getContig(), ci.getStart() + 1, ci.getEnd() + 1, ci.isNegativeStrand(), null);
                        }
                    }

                    if (usedLoci.contains(aci)) {
                        keepGoing = false;
                        break;
                    }

                    String aref = LOOKUPS.get(background).findKmer(aci);

                    CortexRecord cr = GRAPH.findRecord(new CortexKmer(sk));
                    Map<Integer, Set<String>> aks = goForward ? TraversalEngine.getAllNextKmers(cr, !sk.equals(cr.getKmerAsString())) : TraversalEngine.getAllPrevKmers(cr, !sk.equals(cr.getKmerAsString()));
                    Set<String> achi = aks.get(GRAPH.getColorForSampleName(CHILD));

                    if (achi.contains(aref)) {
                        distanceFromNovel++;

                        sk = aref;
                        ci = aci;

                        if (goForward) {
                            vertices.add(sk);
                            loci.add(ci);
                        } else {
                            vertices.add(0, sk);
                            loci.add(0, ci);
                        }
                    } else {
                        if (achi.size() == 1) {
                            distanceFromNovel = 0;
                            sk = achi.iterator().next();
                            Set<Interval> acis = LOOKUPS.get(background).findKmer(sk);
                            ci = acis.size() == 1 ? acis.iterator().next() : null;

                            if ((ci != null && !aci.getContig().equals(ci.getContig())) || acis.size() > 1) {
                                keepGoing = false;
                            } else if (goForward) {
                                vertices.add(sk);
                                loci.add(ci);
                            } else {
                                vertices.add(0, sk);
                                loci.add(0, ci);
                            }
                        } else {
                            keepGoing = false;
                        }

                        onRef = false;
                    }

                    if (!usedLoci.contains(aci)) {
                        usedLoci.add(aci);
                    }
                } while (onRef);
            }
        }

        return new Pair<>(vertices, loci);
    }

    private String recordToString(CortexRecord cr, int childColor, List<Integer> parentColors) {
        List<Object> pieces = new ArrayList<>();
        pieces.add(cr.getKmerAsString());
        pieces.add(cr.getCoverage(childColor));
        parentColors.forEach(c -> pieces.add(cr.getCoverage(c)));
        pieces.add(cr.getEdgesAsString(childColor));
        parentColors.forEach(c -> pieces.add(cr.getEdgesAsString(c)));

        return Joiner.on(' ').join(pieces);
    }
}
