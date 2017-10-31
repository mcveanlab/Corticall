package uk.ac.ox.well.cortexjdk.commands.call.call;

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
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
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
public class Call extends Module {
    @Argument(fullName = "graph", shortName = "g", doc = "Graph")
    public CortexGraph GRAPH;

    @Argument(fullName = "links", shortName = "l", doc = "Links")
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName = "references", shortName = "R", doc = "References")
    public HashMap<String, KmerLookup> REFERENCES;

    @Argument(fullName = "roi", shortName = "r", doc = "ROI")
    public CortexGraph ROI;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CortexKmer, Boolean> seen = new HashMap<>();
        for (CortexRecord cr : ROI) {
            seen.put(cr.getCortexKmer(), false);
        }

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

        Map<CortexKmer, List<CortexVertex>> longWalks = new HashMap<>();

        for (CortexKmer ck : seen.keySet()) {
            if (!seen.get(ck)) {
                List<CortexVertex> l = longWalk(seen, e, ck);

                for (CortexVertex v : l) {
                    if (seen.containsKey(v.getCk())) {
                        if (!longWalks.containsKey(v.getCk()) || longWalks.get(v.getCk()).size() < l.size()) {
                            longWalks.put(v.getCk(), l);
                        }

                        seen.put(v.getCk(), true);
                    }
                }
            }

            pm.update();
        }

        pm = new ProgressMeterFactory()
                .header("Reducing contig set...")
                .message("processed")
                .maxRecord(longWalks.size())
                .make(log);

        Map<String, List<CortexVertex>> longContigs = new HashMap<>();
        for (List<CortexVertex> l : longWalks.values()) {
            String longContig = SequenceUtils.alphanumericallyLowestOrientation(TraversalEngine.toContig(l));
            longContigs.put(longContig, l);

            pm.update();
        }
        log.info("  {} contigs remaining", longContigs.size());

        log.info("Contigs:");
        int contigIndex = 0;
        for (String longContig : longContigs.keySet()) {
            List<CortexVertex> l = longContigs.get(longContig);

            log.info("  {} {} {}", contigIndex, longContig.length(), numNovels(longContigs.get(longContig), seen));

            List<CortexVertex> p = closeBubbles(l, seen);
            log.info("  - novels after bubble closing {}/{}", numNovels(p, seen), numNovels(longContigs.get(longContig), seen));

            if (numNovels(p, seen) > 10) {
                List<List<CortexVertex>> s = breakContigs(p, seen);

                List<SAMRecord> srs = new ArrayList<>();
                Set<String> chrs = new HashSet<>();
                int numPieces = 0;
                int maxLength = 0;
                for (int i = 0; i < s.size(); i++) {
                    List<CortexVertex> q = s.get(i);
                    String contig = TraversalEngine.toContig(q);

                    SAMRecord sr = chooseBestAlignment(contig, 10);
                    srs.add(sr);

                    if (sr != null) {
                        chrs.add(sr.getReferenceName());
                        numPieces++;
                        if (contig.length() > maxLength) {
                            maxLength = contig.length();
                        }
                    }

                }

                if (chrs.size() > 1 && numPieces > 1 && maxLength >= GRAPH.getKmerSize() + 1) {
                    for (int i = 0; i < s.size(); i++) {
                        SAMRecord sr = srs.get(i);

                        log.info("  {} {}", i, sr == null ? "null" : sr.getSAMString().trim());
                    }
                }
            }

            contigIndex++;
        }
    }

    private List<List<CortexVertex>> breakContigs(List<CortexVertex> w, Map<CortexKmer, Boolean> seen) {
        Set<CortexKmer> breakpoints = new HashSet<>();

        boolean inNovelRun = false;
        for (CortexVertex v : w) {
            if (seen.containsKey(v.getCk())) {
                if (!inNovelRun) {
                    breakpoints.add(v.getCk());
                }
                inNovelRun = true;
            } else {
                if (inNovelRun) {
                    breakpoints.add(v.getCk());
                }
                inNovelRun = false;
            }
        }

        List<List<CortexVertex>> s = new ArrayList<>();

        List<CortexVertex> a = new ArrayList<>();
        for (int i = 0; i < w.size(); i++) {
            CortexVertex v = w.get(i);

            if (breakpoints.contains(v.getCk())) {
                s.add(a);
                a = new ArrayList<>();
            }

            a.add(v);
        }

        if (a.size() > 0) {
            s.add(a);
        }

        List<List<CortexVertex>> r = new ArrayList<>();
        for (List<CortexVertex> q : s) {
            if (q.size() > 0 && !seen.containsKey(q.get(0).getCk())) {
                r.add(q);
            }
        }

        return r;
    }

    private int numNovels(List<CortexVertex> w, Map<CortexKmer, Boolean> seen) {
        int numNovels = 0;

        for (CortexVertex v : w) {
            if (seen.containsKey(v.getCk())) {
                numNovels++;
            }
        }

        return numNovels;
    }

    @NotNull
    private List<CortexVertex> longWalk(Map<CortexKmer, Boolean> seen, TraversalEngine e, CortexKmer ck) {
        List<CortexVertex> w = e.walk(ck.getKmerAsString());

        boolean extended;
        do {
            extended = false;
            List<List<CortexVertex>> extFwd = new ArrayList<>();

            Set<CortexVertex> nvs = e.getNextVertices(w.get(w.size() - 1).getBk());
            for (CortexVertex cv : nvs) {
                List<CortexVertex> wn = e.walk(cv.getSk(), true);
                wn.add(0, cv);

                boolean hasNovels = false;

                for (CortexVertex v : wn) {
                    if (seen.containsKey(v.getCk()) && !seen.get(v.getCk())) {
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
                    if (seen.containsKey(v.getCk())) {
                        seen.put(v.getCk(), true);
                    }
                }
            }
        } while (extended);

        do {
            extended = false;
            List<List<CortexVertex>> extRev = new ArrayList<>();

            Set<CortexVertex> pvs = e.getPrevVertices(w.get(0).getBk());
            for (CortexVertex cv : pvs) {
                List<CortexVertex> wp = e.walk(cv.getSk(), false);
                wp.add(cv);

                boolean hasNovels = false;

                for (CortexVertex v : wp) {
                    if (seen.containsKey(v.getCk()) && !seen.get(v.getCk())) {
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
                    if (seen.containsKey(v.getCk())) {
                        seen.put(v.getCk(), true);
                    }
                }
            }
        } while (extended);

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

        @Override
        public String toString() {
            return "LittleBubble{" +
                    "refContig='" + refContig + '\'' +
                    ", altContig='" + altContig + '\'' +
                    ", refPath=" + refPath +
                    ", altPath=" + altPath +
                    ", start=" + start +
                    ", stop=" + stop +
                    '}';
        }
    }

    private List<CortexVertex> closeBubbles(List<CortexVertex> w, Map<CortexKmer, Boolean> seen) {
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
                .joiningColors(GRAPH.getColorForSampleName(ROI.getSampleName(0)))
                .traversalDirection(FORWARD)
                .combinationOperator(OR)
                .stoppingRule(BubbleClosingStopper.class)
                .graph(GRAPH)
                .links(LINKS)
                .references(REFERENCES.values())
                .make();

        Map<Integer, LittleBubble> l = new TreeMap<>();

        for (String parent : REFERENCES.keySet()) {
            e.getConfiguration().setTraversalColor(GRAPH.getColorForSampleName(parent));
            e.getConfiguration().setRecruitmentColors(GRAPH.getColorForSampleName(REFERENCES.get(parent).getSources().iterator().next()));

            for (int i = 0; i < w.size() - 1; i++) {
                CortexVertex vi = w.get(i);

                if (seen.containsKey(vi.getCk())) {
                    List<CortexVertex> roots = new ArrayList<>();
                    List<CortexVertex> sources = new ArrayList<>();

                    int lowerLimit = i - 3 * vi.getCr().getKmerSize() >= 0 ? i - 3 * vi.getCr().getKmerSize() : 0;
                    for (int j = i - 1; j >= lowerLimit; j--) {
                        CortexVertex vj = w.get(j);
                        CortexVertex vk = w.get(j + 1);

                        if (!seen.containsKey(vj.getCk())) {
                            Set<CortexVertex> nvs = e.getNextVertices(new CortexByteKmer(vj.getSk()));

                            for (CortexVertex cv : nvs) {
                                if (!cv.equals(vk)) {
                                    roots.add(vj);
                                    sources.add(cv);
                                }
                            }
                        }
                    }

                    DirectedWeightedPseudograph<CortexVertex, CortexEdge> sinks = new DirectedWeightedPseudograph<>(CortexEdge.class);
                    int distanceFromLastNovel = 0;
                    for (int j = i + 1; j < w.size() && distanceFromLastNovel < 3 * vi.getCr().getKmerSize(); j++) {
                        CortexVertex vj = w.get(j);
                        sinks.addVertex(vj);

                        if (seen.containsKey(vj.getCk())) {
                            distanceFromLastNovel = 0;
                        } else {
                            distanceFromLastNovel++;
                        }
                    }

                    e.getConfiguration().setPreviousTraversal(sinks);

                    for (int q = 0; q < sources.size(); q++) {
                        CortexVertex root = roots.get(q);
                        CortexVertex source = sources.get(q);

                        DirectedWeightedPseudograph<CortexVertex, CortexEdge> b = e.dfs(source.getSk());

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
        }

        List<CortexVertex> wp = new ArrayList<>();
        int usedBubbles = 0;

        for (int i = 0; i < w.size(); i++) {
            if (!l.containsKey(i) || l.get(i).stop < i) {
                wp.add(w.get(i));
            } else {
                LittleBubble lb = l.get(i);

                usedBubbles++;

                for (CortexVertex v : lb.refPath) {
                    wp.add(v);
                }

                i = lb.stop;
            }
        }

        log.info("  - closed {}/{} bubbles", usedBubbles, l.size());

        return wp;
    }

    private SAMRecord chooseBestAlignment(String contig, int mqThreshold) {
        SAMRecord bestRecord = null;
        int bestRecordScore = Integer.MAX_VALUE;

        for (KmerLookup kl : REFERENCES.values()) {
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
