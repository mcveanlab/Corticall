package uk.ac.ox.well.cortexjdk.commands.discover.call;

import com.google.common.base.Joiner;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
import org.jgrapht.GraphPath;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.HTreeMap;
import org.mapdb.Serializer;
import org.mapdb.serializer.SerializerArray;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.mosaic.MosaicAligner;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.DestinationStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ExplorationStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.REVERSE;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalUtils.getAllNextKmers;

public class Call extends Module {
    @Argument(fullName = "graph", shortName = "g", doc = "Graph")
    public CortexGraph GRAPH;

    @Argument(fullName = "links", shortName = "l", doc = "Links", required=false)
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName = "rois", shortName = "r", doc = "Rois")
    public CortexGraph ROIS;

    @Argument(fullName="db", shortName="d", doc="Database")
    public File DB_FILE;

    @Argument(fullName="mother", shortName="m", doc="Mother's sample name")
    public LinkedHashSet<String> MOTHER;

    @Argument(fullName="father", shortName="f", doc="Father's sample name")
    public LinkedHashSet<String> FATHER;

    @Argument(fullName="background", shortName="b", doc="Background", required=false)
    public HashMap<String, IndexedReference> BACKGROUNDS;

    @Argument(fullName="reference", shortName="R", doc="Reference", required=false)
    public HashMap<String, IndexedReference> REFERENCE;

    //@Argument(fullName="reads", shortName="r", doc="Read sequences")
    //public ArrayList<File> READ_FILES;

    @Output
    public PrintStream out;

    @Output(fullName="accOut", shortName="ao", doc="Accounting out")
    public PrintStream aout;

    @Override
    public void execute() {
        log.info("Loading contigs from database...");
        DB dbIn = DBMaker.fileDB(DB_FILE).readOnly().fileMmapEnable().make();

        HTreeMap<Integer, CortexVertex[]> contigsMap = initializeDbContigsMap(dbIn);

        Set<CanonicalKmer> rois = loadRois(ROIS);
        Map<CanonicalKmer, Set<String>> usedRois = new HashMap<>();
        rois.forEach(c -> usedRois.put(c, new TreeSet<>()));

        int sampleColor = GRAPH.getColorForSampleName(ROIS.getSampleName(0));

        Set<Integer> contigIndices = new TreeSet<>(contigsMap.getKeys());
        for (Integer contigIndex : contigIndices) {
            List<CortexVertex> w = getWalk(contigsMap, contigIndex);
            List<CortexVertex> ws = subsetContig(rois, w, 100);

            String contig = TraversalUtils.toContig(ws);

            log.info("contig {}/{}: {} ({})", contigIndex, contigsMap.size(), w.size(), ws.size());

            /*
            List<Pair<Integer, Integer>> regions = getRegions(rois, ws);
            Map<String, String> asmTracks = new HashMap<>();
            for (Set<String> css : Arrays.asList(MOTHER, FATHER)) {
                List<CortexVertex> lv = getParentalContig(sampleColor, ws, regions, css);

                String lc = TraversalUtils.toContig(lv);

                if (!lc.equals(contig)) {
                    asmTracks.put(css.iterator().next(), TraversalUtils.toContig(lv));
                }
            }
            List<Pair<String, String>> aAlignments = ma.align(contig, asmTracks);
            log.info("\n{}", ma);
            */

            MosaicAligner ma = new MosaicAligner();

            List<Pair<String, Interval>> parentIntervals = IntervalCombiner.getIntervals(ws, BACKGROUNDS, 100, 10);
            Map<String, String> parentTracks = loadSequences(parentIntervals, BACKGROUNDS);
            List<Pair<String, String>> pAlignments = ma.align(contig, parentTracks);
            //log.info("\n{}", ma);

            StringBuilder sb = new StringBuilder(StringUtil.repeatCharNTimes(' ', pAlignments.get(0).getSecond().length()));
            int numKmersMarked = 0;
            for (int i = 0; i <= pAlignments.get(0).getSecond().length() - GRAPH.getKmerSize(); i++) {
                int gapsize = GRAPH.getKmerSize() - pAlignments.get(0).getSecond().substring(i, i + GRAPH.getKmerSize()).replaceAll("-", "").length();

                if (i + GRAPH.getKmerSize() + gapsize <= pAlignments.get(0).getSecond().length() - GRAPH.getKmerSize()) {
                    CanonicalKmer ck = new CanonicalKmer(pAlignments.get(0).getSecond().substring(i, i + GRAPH.getKmerSize() + gapsize).replaceAll("-", ""));

                    if (rois.contains(ck)) {
                        numKmersMarked++;
                        for (int j = i; j < i + GRAPH.getKmerSize() + gapsize; j++) {
                            sb.setCharAt(j, 'v');
                        }
                    }
                }
            }

            if (numKmersMarked == 0) {
                sb = new StringBuilder(StringUtil.repeatCharNTimes('v', pAlignments.get(0).getSecond().length()));
            }

            log.info("{} {}", sb.toString(), numKmersMarked);
            for (Pair<String, String> p : pAlignments) {
                log.info("{} {}", p.getSecond(), p.getFirst());
            }

            List<Triple<Integer, String, String>> variants = new ArrayList<>();
            List<Set<CanonicalKmer>> allNovelKmers = new ArrayList<>();

            int variantStart = -1;
            StringBuilder childAllele = new StringBuilder();
            StringBuilder parentAllele = new StringBuilder();
            Set<CanonicalKmer> novelKmers = new HashSet<>();

            String childContig = pAlignments.get(0).getSecond();
            for (int j = 1; j < pAlignments.size(); j++) {
                String parentContig = pAlignments.get(j).getSecond();

                for (int i = 0; i < parentContig.length(); i++) {
                    if (sb.charAt(i) == 'v') {
                        if (childContig.charAt(i) != parentContig.charAt(i) && parentContig.charAt(i) != ' ') {
                            if (variantStart == -1) {
                                variantStart = i;

                                if (i > 0 && (childContig.charAt(i) == '-' || parentContig.charAt(i) == '-')) {
                                    variantStart--;
                                    childAllele.append(childContig.charAt(i-1));
                                    parentAllele.append(parentContig.charAt(i-1));
                                }
                            }

                            if (childContig.charAt(i) != '-') { childAllele.append(childContig.charAt(i)); }
                            if (parentContig.charAt(i) != '-') { parentAllele.append(parentContig.charAt(i)); }
                        } else {
                            if (variantStart >= 0) {
                                variants.add(Triple.of(variantStart, childAllele.toString(), parentAllele.toString()));

                                for (int k = variantStart - GRAPH.getKmerSize(); k < variantStart + childAllele.length() + GRAPH.getKmerSize(); k++) {
                                    if (k >= 0 && k <= childContig.length() - GRAPH.getKmerSize()) {
                                        //log.info("debug: {} {}", k, childContig.length());
                                        String qrs = childContig.substring(k, childContig.length() - 1);
                                        qrs = qrs.replaceAll("-", "");
                                        if (qrs.length() > GRAPH.getKmerSize()) {
                                            CanonicalKmer ck = new CanonicalKmer(qrs.substring(0, GRAPH.getKmerSize()));

                                            if (rois.contains(ck)) {
                                                novelKmers.add(ck);
                                            }
                                        }
                                    }

                                }

                                allNovelKmers.add(novelKmers);

                                variantStart = -1;
                                childAllele = new StringBuilder();
                                parentAllele = new StringBuilder();
                                novelKmers = new HashSet<>();
                            }
                        }
                    }
                }
            }

            int q = 0;
            for (int j = 1; j < pAlignments.size() - 1; j++) {
                q += pAlignments.get(j).getSecond().length();
                variants.add(Triple.of(q, pAlignments.get(j).getFirst(), pAlignments.get(j+1).getFirst()));

                for (int k = q - GRAPH.getKmerSize(); k < q + childAllele.length() + GRAPH.getKmerSize(); k++) {
                    if (k >= 0 && k <= childContig.length() - GRAPH.getKmerSize()) {
                        String qrs = childContig.substring(k, contig.length()).replaceAll("-", "");
                        if (qrs.length() > GRAPH.getKmerSize()) {
                            CanonicalKmer ck = new CanonicalKmer(qrs.substring(0, GRAPH.getKmerSize()));

                            if (rois.contains(ck)) {
                                novelKmers.add(ck);
                            }
                        }
                    }
                }

                allNovelKmers.add(novelKmers);
            }

            if (variantStart >= 0) {
                variants.add(Triple.of(variantStart, childAllele.toString(), parentAllele.toString()));

                for (int k = variantStart - GRAPH.getKmerSize(); k < variantStart + childAllele.length() + GRAPH.getKmerSize(); k++) {
                    if (k >= 0 && k <= childContig.length() - GRAPH.getKmerSize()) {
                        String qrs = childContig.substring(k, contig.length()).replaceAll("-", "");
                        if (qrs.length() > GRAPH.getKmerSize()) {
                            CanonicalKmer ck = new CanonicalKmer(qrs.substring(0, GRAPH.getKmerSize()));

                            if (rois.contains(ck)) {
                                novelKmers.add(ck);
                            }
                        }
                    }
                }

                allNovelKmers.add(novelKmers);
            }

            for (int i = 0; i < variants.size(); i++) {
                Triple<Integer, String, String> variant = variants.get(i);
                log.info("  {} {}", variant, allNovelKmers.get(i).size());

                for (CanonicalKmer ck : allNovelKmers.get(i)) {
                    usedRois.get(ck).add(contigIndex + ":" + i);
                }

                out.println(Joiner.on("\t").join(contigIndex, i, variant.getLeft(), variant.getMiddle(), variant.getRight(), Joiner.on(",").join(allNovelKmers.get(i))));
            }

            for (CanonicalKmer ck : usedRois.keySet()) {
                if (usedRois.get(ck).size() == 0) {
                    aout.println(ck + "\tNA");
                } else {
                    aout.println(ck + "\t" + Joiner.on(",").join(usedRois.get(ck)));
                }
            }

            /*
            for (String refid : REFERENCE.keySet()) {
                Map<String, IndexedReference> ref = new HashMap<>();
                ref.put(refid, REFERENCE.get(refid));

                List<Pair<String, Interval>> refIntervals = IntervalCombiner.getIntervals(ws, ref, 100, 10);
                Map<String, String> refTracks = loadSequences(refIntervals, ref);
                List<Pair<String, String>> rAlignments = ma.align(contig, refTracks);
                log.info("\n{}", ma);
            }
            */
        }
    }

    private List<CortexVertex> getParentalContig(int sampleColor, List<CortexVertex> ws, List<Pair<Integer, Integer>> regions, Set<String> css) {
        List<Triple<Integer, Integer, List<CortexVertex>>> replacements = new ArrayList<>();

        for (Pair<Integer, Integer> r : regions) {
            Set<String> left = getLeft(sampleColor, ws, r, css);
            Set<String> right = getRight(sampleColor, ws, r, css);

            TraversalEngine ef = new TraversalEngineFactory()
                    .traversalColors(GRAPH.getColorsForSampleNames(css))
                    .combinationOperator(OR)
                    .traversalDirection(FORWARD)
                    .graph(GRAPH)
                    .links(LINKS)
                    .stoppingRule(DestinationStopper.class)
                    .make();

            TraversalEngine er = new TraversalEngineFactory()
                    .traversalColors(GRAPH.getColorsForSampleNames(css))
                    .combinationOperator(OR)
                    .traversalDirection(REVERSE)
                    .graph(GRAPH)
                    .links(LINKS)
                    .stoppingRule(DestinationStopper.class)
                    .make();

            DirectedWeightedPseudograph<CortexVertex, CortexEdge> gf = ef.dfs(left, right);
            DirectedWeightedPseudograph<CortexVertex, CortexEdge> gr = er.dfs(right, left);

            for (DirectedWeightedPseudograph<CortexVertex, CortexEdge> g : Arrays.asList(gf, gr)) {
                if (g != null && g.vertexSet().size() > 0) {
                    Triple<Integer, Integer, List<CortexVertex>> replacement = getReplacement(ws, left, right, g);

                    if (replacement != null) {
                        replacements.add(replacement);
                        break;
                    }
                }
            }
        }

        List<CortexVertex> qss = new ArrayList<>();

        Iterator<Triple<Integer, Integer, List<CortexVertex>>> rit = replacements.iterator();
        Triple<Integer, Integer, List<CortexVertex>> currentReplacement = rit.hasNext() ? rit.next() : null;
        for (int i = 0; i < ws.size(); i++) {
            if (currentReplacement == null || i < currentReplacement.getLeft()) {
                qss.add(ws.get(i));
            } else if (i == currentReplacement.getLeft()) {
                qss.addAll(currentReplacement.getRight());
                i = currentReplacement.getMiddle();

                currentReplacement = rit.hasNext() ? rit.next() : null;
            }
        }

        return qss;
    }

    private Triple<Integer, Integer, List<CortexVertex>> getReplacement(List<CortexVertex> ws, Set<String> left, Set<String> right, DirectedWeightedPseudograph<CortexVertex, CortexEdge> g) {
        PathFinder pf = new PathFinder(g, g.edgeSet().iterator().next().getColor());

        for (String before : left) {
            for (String after : right) {
                List<GraphPath<CortexVertex, CortexEdge>> gps = pf.getPaths(TraversalUtils.findVertex(g, before), TraversalUtils.findVertex(g, after));

                for (GraphPath<CortexVertex, CortexEdge> gp : gps) {
                    int start = 0, end = ws.size() - 1;
                    for (int i = 0; i < ws.size(); i++) {
                        if (ws.get(i).getKmerAsString().equals(before)) { start = i; }
                        if (ws.get(i).getKmerAsString().equals(after)) { end = i; }
                    }

                    return Triple.of(start, end, gp.getVertexList());
                }
            }
        }

        return null;
    }

    private Set<String> getRight(int sampleColor, List<CortexVertex> ws, Pair<Integer, Integer> r, Set<String> css) {
        Set<String> rightKmers = new LinkedHashSet<>();

        int right = 0;
        for (right = r.getSecond() + 1; right < ws.size() - 2 && right <= r.getSecond() + 50; right++) {
            Map<Integer, Set<CortexByteKmer>> cbks = TraversalUtils.getAllPrevKmers(ws.get(right).getCortexRecord(), !ws.get(right).getKmerAsString().equals(ws.get(right).getCanonicalKmer().getKmerAsString()));

            for (int c : GRAPH.getColorsForSampleNames(css)) {
                if (cbks.get(c).size() > 0 && !cbks.get(c).equals(cbks.get(sampleColor))) {
                    rightKmers.add(ws.get(right).getKmerAsString());
                }
            }
        }

        return rightKmers;
    }

    private Set<String> getLeft(int sampleColor, List<CortexVertex> ws, Pair<Integer, Integer> r, Set<String> css) {
        Set<String> leftKmers = new LinkedHashSet<>();

        int left = 0;
        for (left = r.getFirst() - 1; left > 0 && left >= r.getFirst() - 50; left--) {
            Map<Integer, Set<CortexByteKmer>> cbks = TraversalUtils.getAllNextKmers(ws.get(left).getCortexRecord(), !ws.get(left).getKmerAsString().equals(ws.get(left).getCanonicalKmer().getKmerAsString()));

            for (int c : GRAPH.getColorsForSampleNames(css)) {
                if (cbks.get(c).size() > 0 && !cbks.get(c).equals(cbks.get(sampleColor))) {
                    leftKmers.add(ws.get(left).getKmerAsString());
                }
            }
        }

        return leftKmers;
    }

    @NotNull
    private Map<String, String> loadSequences(List<Pair<String, Interval>> bfParents, Map<String, IndexedReference> refs) {
        Map<String, String> targets = new LinkedHashMap<>();
        for (Pair<String, Interval> l : bfParents) {
            Interval it = l.getSecond();

            String seq = refs.get(l.getFirst()).find(it);
            String targetName = l.getFirst() + ":" + it.getContig() + ":" + it.getStart() + "-" + it.getEnd() + ":" + (it.isPositiveStrand() ? "+" : "-");
            targets.put(targetName, seq);
        }
        return targets;
    }

    private List<Triple<Integer, String, String>> call(String[] al) {
        List<Triple<Integer, String, String>> calls = new ArrayList<>();

        StringBuilder cAllele = new StringBuilder();
        StringBuilder pAllele = new StringBuilder();

        for (int i = 0; i < al[0].length(); i++) {
            if (al[0].charAt(i) == al[1].charAt(i)) {
                if (cAllele.length() > 0 || pAllele.length() > 0) {
                    String cStr = cAllele.toString().replaceAll("-", "");
                    String pStr = pAllele.toString().replaceAll("-", "");

                    calls.add(Triple.of(i - cAllele.length(), cStr, pStr));

                    cAllele = new StringBuilder();
                    pAllele = new StringBuilder();
                }
            } else {
                cAllele.append(al[0].charAt(i));
                pAllele.append(al[1].charAt(i));
            }
        }

        if (cAllele.length() > 0 || pAllele.length() > 0) {
            String cStr = cAllele.toString().replaceAll("-", "");
            String pStr = pAllele.toString().replaceAll("-", "");

            calls.add(Triple.of(al[0].length() - cAllele.length(), cStr, pStr));
        }

        return calls;
    }

    private List<CortexVertex> subsetContig(Set<CanonicalKmer> rois, List<CortexVertex> w, int window) {
        List<Pair<Integer, Integer>> regions = getRegions(rois, w);

        int subcontigStart = regions.get(0).getFirst() - window;
        if (subcontigStart < 0) { subcontigStart = 0; }

        int subcontigStop  = regions.get(regions.size() - 1).getSecond() + window;
        if (subcontigStop >= w.size()) { subcontigStop = w.size() - 1; }

        List<CortexVertex> ws = new ArrayList<>();
        for (int i = subcontigStart; i <= subcontigStop; i++) {
            ws.add(w.get(i));
        }

        return ws;
    }

    @NotNull
    private List<Pair<Integer, Integer>> getRegions(Set<CanonicalKmer> rois, List<CortexVertex> cvs) {
        List<Pair<Integer, Integer>> regions = new ArrayList<>();

        int regionStart = -1;
        int regionStop = 0;
        for (int i = 0; i < cvs.size(); i++) {
            CanonicalKmer currentKmer = cvs.get(i).getCanonicalKmer();

            if (rois.contains(currentKmer)) {
                if (regionStart == -1) { regionStart = i; }
                regionStop = i;
            } else {
                if (regionStart > -1) {
                    regions.add(new Pair<>(regionStart, regionStop));

                    regionStart = -1;
                    regionStop = 0;
                }
            }
        }

        if (regionStart > -1) {
            regions.add(new Pair<>(regionStart, regionStop));
        }

        return regions;
    }

    private List<CortexVertex> getWalk(HTreeMap<Integer, CortexVertex[]> contigsMap, int contigIndex) {
        List<CortexVertex> cvs = new ArrayList<>();

        Object[] oa = contigsMap.get(contigIndex);

        for (Object o : oa) {
            cvs.add((CortexVertex) o);
        }

        return cvs;
    }

    /*
    private HTreeMap<Integer, FastqRecord[]> initializeDbReadsMap(DB dbIn, boolean firstEnd) {
        HTreeMap<Integer, FastqRecord[]> dbReads = dbIn
                .hashMap("readsEnd" + (firstEnd ? "1" : "2"))
                .keySerializer(Serializer.INTEGER)
                .valueSerializer(new SerializerArray(Serializer.JAVA))
                .counterEnable()
                .open();

        return dbReads;
    }
    */

    @NotNull
    private HTreeMap<Integer, CortexVertex[]> initializeDbContigsMap(DB dbIn) {
        return (HTreeMap<Integer, CortexVertex[]>) dbIn
                .hashMap("contigs")
                .keySerializer(Serializer.INTEGER)
                .valueSerializer(new SerializerArray(Serializer.JAVA))
                .counterEnable()
                .open();
    }

    private Set<CanonicalKmer> loadRois(CortexGraph rois) {
        Set<CanonicalKmer> roiSet = new HashSet<>();

        for (CortexRecord rr : rois) {
            roiSet.add(rr.getCanonicalKmer());
        }

        return roiSet;
    }
}
