package uk.ac.ox.well.cortexjdk.commands.call.eval;

import com.google.common.base.Joiner;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;
import org.jetbrains.annotations.NotNull;
import org.jgrapht.GraphPath;
import org.jgrapht.alg.shortestpath.DijkstraShortestPath;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.BubbleClosingStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.NovelContinuationStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;

/**
 * Created by kiran on 30/08/2017.
 */
public class CompareContigLengths extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Links")
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName="references", shortName="R", doc="References")
    public HashMap<String, IndexedReference> REFERENCES;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CanonicalKmer, Boolean> seen = new HashMap<>();
        for (CortexRecord cr : ROI) {
            seen.put(cr.getCanonicalKmer(), false);
        }

        TraversalEngine eo = new TraversalEngineFactory()
                .traversalColor(GRAPH.getColorForSampleName(ROI.getSampleName(0)))
                .joiningColors(GRAPH.getColorsForSampleNames(REFERENCES.keySet()))
                .combinationOperator(OR)
                .stoppingRule(NovelContinuationStopper.class)
                .graph(GRAPH)
                .references(REFERENCES.values())
                .rois(ROI)
                .make();

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(GRAPH.getColorForSampleName(ROI.getSampleName(0)))
                .joiningColors(GRAPH.getColorsForSampleNames(REFERENCES.keySet()))
                .combinationOperator(OR)
                .stoppingRule(NovelContinuationStopper.class)
                .graph(GRAPH)
                .links(LINKS)
                .references(REFERENCES.values())
                .rois(ROI)
                .make();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing novel kmers")
                .message("processed")
                .maxRecord(seen.size())
                .make(log);

        Map<CanonicalKmer, List<CortexVertex>> shortWalks = new HashMap<>();
        Map<CanonicalKmer, List<CortexVertex>> longWalks = new HashMap<>();

        out.println(Joiner.on("\t").join("ck", "wsize", "lsize", "w", "l"));
        for (CanonicalKmer ck : seen.keySet()) {
            List<CortexVertex> w = eo.walk(ck.getKmerAsString());
            List<CortexVertex> l = longWalk(seen, e, ck);

            for (CortexVertex v : w) {
                if (seen.containsKey(v.getCanonicalKmer())) {
                    if (!shortWalks.containsKey(v.getCanonicalKmer()) || shortWalks.get(v.getCanonicalKmer()).size() < w.size()) {
                        shortWalks.put(v.getCanonicalKmer(), w);
                    }
                }
            }

            for (CortexVertex v : l) {
                if (!longWalks.containsKey(v.getCanonicalKmer()) || longWalks.get(v.getCanonicalKmer()).size() < l.size()) {
                    longWalks.put(v.getCanonicalKmer(), l);
                }
            }

            out.println(Joiner.on("\t").join(ck, w.size(), l.size(), TraversalEngine.toContig(w), TraversalEngine.toContig(l)));
            log.info("  {} {} {}", ck, w.size(), l.size());

            pm.update();
        }

        /*
        Set<String> shortContigs = new HashSet<>();
        Set<String> longContigs = new HashSet<>();

        for (List<CortexVertex> w : shortWalks.values()) {
            String shortContig = SequenceUtils.alphanumericallyLowestOrientation(TraversalEngine.toContig(w));
            shortContigs.add(shortContig);
        }

        for (List<CortexVertex> l : longWalks.values()) {
            String longContig = SequenceUtils.alphanumericallyLowestOrientation(TraversalEngine.toContig(l));
            longContigs.add(longContig);
        }

        log.info("{} {}", shortContigs.size(), longContigs.size());
        */

        /*
        for (CortexKmer ck : seen.keySet()) {
            if (!seen.get(ck)) {
                List<CortexVertex> w = eo.walk(ck.getKmerAsString());
                List<CortexVertex> l = longWalk(seen, e, ck);

                //log.info("    orig: {}", TraversalEngine.toContig(w));
                //log.info("    open: {}", TraversalEngine.toContig(l));

                aout.println(w.size() + "\t" + l.size());

                log.info("  short walk: {}, long walk: {}, num novels: {}", w.size(), l.size(), numNovels(l, seen));

                List<List<CortexVertex>> contigsWithClosedBubbles = null;
                int numNovelsRemaining = Integer.MAX_VALUE;

                for (String parent : REFERENCES.keySet()) {
                    List<CortexVertex> p = closeBubbles(l, parent, seen);

                    //log.info("  closed: {}", TraversalEngine.toContig(p));

                    //log.info("  num novels: {}", numNovels(p, seen));

                    List<List<CortexVertex>> s = breakContigs(p, parent, seen);

                    //log.info("{} {} {} {} {} {} {} {} {}", ck, parent, w.size(), l.size(), p.size(), numNovels(w, seen), numNovels(l, seen), numNovels(p, seen), s.size());

                    if (numNovels(p, seen) < numNovelsRemaining) {
                        contigsWithClosedBubbles = s;
                        numNovelsRemaining = numNovels(p, seen);
                    }
                }


                if (contigsWithClosedBubbles != null && numNovelsRemaining >= 10) {
                    List<SAMRecord> srs = new ArrayList<>();
                    int numGoodAlignments = 0;

                    for (int i = 0; i < contigsWithClosedBubbles.size(); i++) {
                        List<CortexVertex> p = contigsWithClosedBubbles.get(i);
                        String contig = TraversalEngine.toContig(p);

                        SAMRecord sr = chooseBestAlignment(contig, MQ_THRESHOLD);
                        srs.add(sr);

                        if (sr != null) {
                            numGoodAlignments++;
                        }
                    }

                    if (numGoodAlignments > 1) {
                        for (int i = 0; i < srs.size(); i++) {
                            SAMRecord sr = srs.get(i);
                            //log.info("{}", sr.getSAMString().trim());

                            if (sr == null) {
                                out.println(Joiner.on("\t").join(
                                        ck,
                                        i,
                                        contigsWithClosedBubbles.get(i).size(),
                                        ".",
                                        ".",
                                        ".",
                                        ".",
                                        ".",
                                        "."
                                        )
                                );
                            } else {
                                out.println(Joiner.on("\t").join(
                                        ck,
                                        i,
                                        contigsWithClosedBubbles.get(i).size(),
                                        sr.getReferenceName(),
                                        sr.getAlignmentStart(),
                                        sr.getAlignmentEnd(),
                                        sr.getReadNegativeStrandFlag() ? "-" : "+",
                                        sr.getIntegerAttribute("NM"),
                                        sr.getCigarString()
                                        )
                                );
                            }
                        }
                    }
                }

                for (CortexVertex v : l) {
                    if (seen.containsKey(v.getCanonicalKmer())) {
                        seen.put(v.getCanonicalKmer(), true);
                        pm.update();
                    }
                }
            }
        }
        */
    }

    private List<List<CortexVertex>> breakContigs(List<CortexVertex> w, String parent, Map<CanonicalKmer, Boolean> seen) {
        Set<CanonicalKmer> breakpoints = new HashSet<>();
        List<Set<Interval>> intervals = new ArrayList<>();

        boolean inNovelRun = false;
        for (CortexVertex v : w) {
            if (seen.containsKey(v.getCanonicalKmer())) {
                if (!inNovelRun) {
                    breakpoints.add(v.getCanonicalKmer());
                }
                inNovelRun = true;
            } else {
                inNovelRun = false;
            }

            intervals.add(REFERENCES.get(parent).find(v.getKmerAsString()));
        }

        for (int i = 0; i < intervals.size(); i++) {
            if (intervals.get(i).size() == 2) {
                Interval leftInterval = null;
                int leftIndex = -1;
                for (int j = i - 1; j >= 0; j--) {
                    if (intervals.get(j).size() == 1) {
                        leftInterval = intervals.get(j).iterator().next();
                        leftIndex = j;
                        break;
                    }
                }

                Interval rightInterval = null;
                int rightIndex = -1;
                int lastMultiAlignmentIndex = -1;
                for (int j = i; j < intervals.size(); j++) {
                    if (intervals.get(j).size() > 1) {
                        lastMultiAlignmentIndex = j;
                    } else if (intervals.get(j).size() == 1) {
                        rightInterval = intervals.get(j).iterator().next();
                        rightIndex = j;
                        break;
                    }
                }

                if (leftInterval != null && rightInterval != null && !leftInterval.intersects(rightInterval)) {
                    boolean leftFound = false;
                    for (Interval it : intervals.get(i)) {
                        if (leftInterval.intersects(it)) {
                            leftFound = true;
                        }
                    }

                    boolean rightFound = false;
                    for (Interval it : intervals.get(lastMultiAlignmentIndex)) {
                        if (rightInterval.intersects(it)) {
                            rightFound = true;
                        }
                    }

                    if (leftFound && rightFound) {
                        breakpoints.add(w.get(i).getCanonicalKmer());

                        i = rightIndex;
                    }
                }
            }
        }

        List<List<CortexVertex>> s = new ArrayList<>();

        List<CortexVertex> a = new ArrayList<>();
        for (CortexVertex v : w) {
            a.add(v);

            if (breakpoints.contains(v.getCanonicalKmer())) {
                s.add(a);
                a = new ArrayList<>();
            }
        }

        if (a.size() > 0) {
            s.add(a);
        }

        return s;
    }

    private int numNovels(List<CortexVertex> w, Map<CanonicalKmer, Boolean> seen) {
        int numNovels = 0;

        for (CortexVertex v : w) {
            if (seen.containsKey(v.getCanonicalKmer())) {
                numNovels++;
            }
        }

        return numNovels;
    }

    @NotNull
    private List<CortexVertex> longWalk(Map<CanonicalKmer, Boolean> seen, TraversalEngine e, CanonicalKmer ck) {
        List<CortexVertex> w = e.walk(ck.getKmerAsString());

        //int wOldSize = w.size();

        boolean extended = false;
        do {
            extended = false;
            List<List<CortexVertex>> extFwd = new ArrayList<>();

            Set<CortexVertex> nvs = e.getNextVertices(w.get(w.size() - 1).getKmerAsByteKmer());
            for (CortexVertex cv : nvs) {
                List<CortexVertex> wn = e.walk(cv.getKmerAsString(), true);
                wn.add(0, cv);

                boolean hasNovels = false;

                for (CortexVertex v : wn) {
                    if (seen.containsKey(v.getCanonicalKmer()) && !seen.get(v.getCanonicalKmer())) {
                        hasNovels = true;
                        break;
                    }
                }

                if (hasNovels) {
                    extFwd.add(wn);
                }
            }

            if (extFwd.size() == 1) {
                w.addAll(extFwd.get(0));
                extended = true;

                for (CortexVertex v : extFwd.get(0)) {
                    if (seen.containsKey(v.getCanonicalKmer())) {
                        seen.put(v.getCanonicalKmer(), true);
                    }
                }
            }
        } while (extended);

        do {
            extended = false;
            List<List<CortexVertex>> extRev = new ArrayList<>();

            Set<CortexVertex> pvs = e.getPrevVertices(w.get(0).getKmerAsByteKmer());
            for (CortexVertex cv : pvs) {
                List<CortexVertex> wp = e.walk(cv.getKmerAsString(), false);
                wp.add(cv);

                boolean hasNovels = false;

                for (CortexVertex v : wp) {
                    if (seen.containsKey(v.getCanonicalKmer()) && !seen.get(v.getCanonicalKmer())) {
                        hasNovels = true;
                        break;
                    }
                }

                if (hasNovels) {
                    extRev.add(wp);
                }
            }

            if (extRev.size() == 1) {
                w.addAll(0, extRev.get(0));
                extended = true;

                for (CortexVertex v : extRev.get(0)) {
                    if (seen.containsKey(v.getCanonicalKmer())) {
                        seen.put(v.getCanonicalKmer(), true);
                    }
                }
            }
        } while (extended);

        //int wNewSize = w.size();

        return w;
    }

    private static String[] contigsToAlleles(String s0, String s1) {
        int s0start = 0, s0end = s0.length();
        int s1start = 0, s1end = s1.length();

        for (int i = 0, j = 0; i < s0.length() && j < s1.length(); i++, j++) {
            if (s0.charAt(i) != s1.charAt(j)) {
                s0start = i;
                s1start = j;
                break;
            }
        }

        for (int i = s0.length() - 1, j = s1.length() - 1; i >= 0 && j >= 0; i--, j--) {
            if (s0.charAt(i) != s1.charAt(j) || i == s0start - 1 || j == s1start - 1) {
                s0end = i + 1;
                s1end = j + 1;
                break;
            }
        }

        String[] pieces = new String[4];
        pieces[0] = s0.substring(0, s0start);
        pieces[1] = s0.substring(s0start, s0end);
        pieces[2] = s1.substring(s1start, s1end);
        pieces[3] = s0.substring(s0end, s0.length());

        return pieces;
    }

    private class LittleBubble implements Comparable<LittleBubble> {
        public String refContig;
        public String altContig;
        public List<CortexVertex> refPath;
        public List<CortexVertex> altPath;
        public Integer start;
        public Integer stop;

        @Override
        public int compareTo(@NotNull LittleBubble o) {
            return start.compareTo(o.start);
        }
    }

    private List<CortexVertex> closeBubbles(List<CortexVertex> w, String parent, Map<CanonicalKmer, Boolean> seen) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = new DirectedWeightedPseudograph<>(CortexEdge.class);
        Map<CortexVertex, Integer> indices = new HashMap<>();
        indices.put(w.get(0), 0);

        for (int i = 1; i < w.size(); i++) {
            CortexVertex v0 = w.get(i - 1);
            CortexVertex v1 = w.get(i);

            g.addVertex(v0);
            g.addVertex(v1);
            g.addEdge(v0, v1, new CortexEdge());

            indices.put(w.get(i), i);
        }

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(GRAPH.getColorForSampleName(parent))
                .joiningColors(GRAPH.getColorForSampleName(ROI.getSampleName(0)))
                .recruitmentColors(GRAPH.getColorForSampleName(REFERENCES.get(parent).getSources().iterator().next()))
                .traversalDirection(FORWARD)
                .combinationOperator(OR)
                .stoppingRule(BubbleClosingStopper.class)
                //.previousTraversal(g)
                .graph(GRAPH)
                .links(LINKS)
                .references(REFERENCES.values())
                .make();

        Map<Integer, LittleBubble> l = new TreeMap<>();

        for (int i = 0; i < w.size() - 1; i++) {
            CortexVertex vi = w.get(i);

            if (seen.containsKey(vi.getCanonicalKmer())) {
                List<CortexVertex> roots = new ArrayList<>();
                List<CortexVertex> sources = new ArrayList<>();

                int lowerLimit = i - 3*vi.getCortexRecord().getKmerSize() >= 0 ? i - 3*vi.getCortexRecord().getKmerSize() : 0;
                for (int j = i - 1; j >= lowerLimit; j--) {
                    CortexVertex vj = w.get(j);
                    CortexVertex vk = w.get(j+1);

                    if (!seen.containsKey(vj.getCanonicalKmer())) {
                        Set<CortexVertex> nvs = e.getNextVertices(new CortexByteKmer(vj.getKmerAsString()));

                        for (CortexVertex cv : nvs) {
                            if (!cv.equals(vk)) {
                                roots.add(vj);
                                sources.add(cv);
                            }
                        }
                    }
                }

                DirectedWeightedPseudograph<CortexVertex, CortexEdge> sinks = new DirectedWeightedPseudograph<>(CortexEdge.class);
                for (int j = i + 1; j < w.size(); j++) {
                    CortexVertex vj = w.get(j);
                    sinks.addVertex(vj);
                }
                e.getConfiguration().setPreviousTraversal(sinks);

                for (int q = 0; q < sources.size(); q++) {
                    CortexVertex root = roots.get(q);
                    CortexVertex source = sources.get(q);

                    DirectedWeightedPseudograph<CortexVertex, CortexEdge> b = e.dfs(source.getKmerAsString());

                    if (b != null) {
                        CortexVertex sink = null;
                        int sinkIndex = -1;
                        for (CortexVertex v : b.vertexSet()) {
                            if (indices.containsKey(v) && indices.get(v) > sinkIndex) {
                                sink = v;
                                sinkIndex = indices.get(sink);
                            }
                        }

                        GraphPath<CortexVertex, CortexEdge> gp = DijkstraShortestPath.findPathBetween(b, source, sink);

                        List<CortexVertex> refPath = new ArrayList<>();
                        refPath.add(root);

                        for (CortexVertex v : gp.getVertexList()) {
                            refPath.add(v);
                        }

                        List<CortexVertex> altPath = new ArrayList<>();

                        for (int j = indices.get(root); j <= indices.get(sink); j++) {
                            altPath.add(w.get(j));
                        }

                        String refContig = TraversalEngine.toContig(refPath);
                        String altContig = TraversalEngine.toContig(altPath);

                        LittleBubble lb = new LittleBubble();
                        lb.refContig = refContig;
                        lb.altContig = altContig;
                        lb.refPath = refPath;
                        lb.altPath = altPath;
                        lb.start = indices.get(root);
                        lb.stop = indices.get(sink);

                        l.put(lb.start, lb);

                        i = lb.stop - 1;
                    }
                }
            }
        }

        for (LittleBubble lb : l.values()) {
            //log.info("lb: {}", lb);
            List<SAMRecord> srs = REFERENCES.get(parent).getAligner().align(lb.refContig);
            for (SAMRecord sr : srs) {
                if (sr.getMappingQuality() > 0) {
                    //log.info("  {}", sr.getSAMString().trim());
                }
            }
        }

        List<CortexVertex> wp = new ArrayList<>();

        for (int i = 0; i < w.size(); i++) {
            if (!l.containsKey(i)) {
                wp.add(w.get(i));
            } else {
                LittleBubble lb = l.get(i);
                for (CortexVertex v : lb.refPath) {
                    wp.add(v);
                }

                i = lb.stop;
            }
        }

        //log.info("  closed {} bubbles", l.size());

        return wp;
    }

    private SAMRecord chooseBestAlignment(String contig, int mqThreshold) {
        SAMRecord bestRecord = null;
        int bestRecordScore = Integer.MAX_VALUE;

        for (IndexedReference kl : REFERENCES.values()) {
            List<SAMRecord> recs = kl.getAligner().align(contig);
            List<SAMRecord> filteredRecs = new ArrayList<>();

            for (SAMRecord rec : recs) {
                if (rec.getMappingQuality() >= mqThreshold) {
                    filteredRecs.add(rec);
                }
            }

            if (filteredRecs.size() == 1) {
                SAMRecord filteredRec = filteredRecs.get(0);
                int filteredRecScore = scoreAlignment(filteredRec);

                if (filteredRecScore < bestRecordScore) {
                    bestRecord = filteredRec;
                    bestRecordScore = filteredRecScore;
                }
            }
        }

        return bestRecord;
    }

    private int scoreAlignment(SAMRecord record) {
        int basesChanged = record.getIntegerAttribute("NM");

        for (CigarElement ce : record.getCigar()) {
            if (ce.getOperator().equals(CigarOperator.SOFT_CLIP) ||
                ce.getOperator().equals(CigarOperator.INSERTION) ||
                ce.getOperator().equals(CigarOperator.DELETION)) {
                basesChanged += ce.getLength();
            }
        }

        return basesChanged;
    }
}
