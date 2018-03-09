package uk.ac.ox.well.cortexjdk.commands.discover.call;

import com.google.common.base.Joiner;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
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
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.alignment.sw.SmithWaterman;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.commands.discover.call.BackgroundInference.BackgroundGroup.CHILD;
import static uk.ac.ox.well.cortexjdk.commands.discover.call.BackgroundInference.BackgroundGroup.PARENT;
import static uk.ac.ox.well.cortexjdk.commands.discover.call.BackgroundInference.BackgroundGroup.REF;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.REVERSE;

public class Annotate extends Module {
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
    public File out;

    @Output(fullName="variantsOut", shortName="vo", doc="Variants out")
    public PrintStream vout;

    @Override
    public void execute() {
        log.info("Loading contigs from database...");
        DB dbIn = DBMaker.fileDB(DB_FILE).readOnly().fileMmapEnable().make();

        HTreeMap<Integer, CortexVertex[]> contigsMap = initializeDbContigsMap(dbIn);
        //HTreeMap<Integer, FastqRecord[]> readsMapEnd1 = initializeDbReadsMap(dbIn, true);
        //HTreeMap<Integer, FastqRecord[]> readsMapEnd2 = initializeDbReadsMap(dbIn, false);

        String sample = ROIS.getSampleName(0);
        Set<CanonicalKmer> rois = loadRois(ROIS);

        Set<Integer> contigIndices = new TreeSet<>(contigsMap.getKeys());
        for (Integer contigIndex : contigIndices) {
            List<CortexVertex> w = getWalk(contigsMap, contigIndex);
            List<CortexVertex> ws = subsetContig(rois, w, 100);

            log.info("contig {}/{}: {} ({})", contigIndex, contigsMap.size(), w.size(), ws.size());

            BackgroundInference bf = new BackgroundInference();
            for (int i = 0; i < ws.size(); i++) {
                CortexVertex v = ws.get(i);
                boolean isNovel = rois.contains(v.getCanonicalKmer());

                bf.add(i, v, sample, new TreeSet<>(), true, isNovel, CHILD);

                if (MOTHER != null && FATHER != null) {
                    for (Set<String> bgNames : Arrays.asList(MOTHER, FATHER)) {
                        boolean isPresent = false;
                        Set<Interval> locations = new TreeSet<>();

                        for (String bgName : bgNames) {
                            isPresent |= v.getCortexRecord().getCoverage(GRAPH.getColorForSampleName(bgName)) > 0;

                            if (BACKGROUNDS != null && BACKGROUNDS.containsKey(bgName)) {
                                locations.addAll(BACKGROUNDS.get(bgName).find(v.getKmerAsString()));
                            }
                        }

                        bf.add(i, v, bgNames.iterator().next(), locations, isPresent, isNovel, PARENT);
                    }
                }

                if (REFERENCE != null) {
                    for (String refName : REFERENCE.keySet()) {
                        boolean isPresent = v.getCortexRecord().getCoverage(GRAPH.getColorForSampleName(refName)) > 0;
                        Set<Interval> locations = new TreeSet<>(REFERENCE.get(refName).find(v.getKmerAsString()));

                        bf.add(i, v, refName, locations, isPresent, isNovel, REF);
                    }
                }
            }

            for (BackgroundInference.BackgroundGroup bg : Arrays.asList(PARENT)) {
                Triple<List<String>, double[][][], List<String>> f = bf.prepare(bg);
                List<String> stateList = f.getLeft();
                double[][][] a = f.getMiddle();
                List<String> stateSeq = f.getRight();

                log.info("  {} {}", stateList.size(), Joiner.on(",").join(stateList));

                //printContigMatrix(contigIndex, rois, ws, a, stateList, stateSeq, bg.toString());

                List<Pair<Integer, Integer>> regions = getRegions(rois, ws);

                List<Integer> acolors = new ArrayList<>();
                if (bg.equals(PARENT)) {
                    acolors.addAll(GRAPH.getColorsForSampleNames(MOTHER));
                    acolors.addAll(GRAPH.getColorsForSampleNames(FATHER));
                } else {
                    acolors.addAll(GRAPH.getColorsForSampleNames(REFERENCE.keySet()));
                }

                for (Pair<Integer, Integer> region : regions) {
                    String childContig = "";

                    int offset = 0;
                    int finalColor = -1;
                    String contig = "";
                    int finalDovetail = 0;
                    String finalOverlap = "";

                    for (int c : acolors) {
                        TraversalEngine ef = new TraversalEngineFactory()
                                .traversalColors(c)
                                .traversalDirection(FORWARD)
                                .combinationOperator(OR)
                                .stoppingRule(ContigStopper.class)
                                .graph(GRAPH)
                                .links(LINKS)
                                .make();

                        TraversalEngine er = new TraversalEngineFactory()
                                .traversalColors(c)
                                .traversalDirection(REVERSE)
                                .combinationOperator(OR)
                                .stoppingRule(ContigStopper.class)
                                .graph(GRAPH)
                                .links(LINKS)
                                .make();

                        for (int expand = 1; expand < 10 && contig.isEmpty(); expand++) {
                            int left = region.getFirst() - expand;
                            if (left < 0) { left = 0; }

                            int right = region.getSecond() + expand;
                            if (right >= ws.size()) { right = ws.size() - 1; }

                            String source = ws.get(left).getKmerAsString();
                            String sink = ws.get(right).getKmerAsString();

                            List<CortexVertex> wss = new ArrayList<>();
                            for (int j = left; j <= right; j++) {
                                wss.add(ws.get(j));
                            }

                            childContig = TraversalUtils.toContig(wss);

                            DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = null;
                            if (g == null) { g = ef.dfs(source, sink); }
                            if (g == null) { g = er.dfs(sink, source); }

                            if (g != null) {
                                CortexVertex cvSource = TraversalUtils.findVertex(g, source);
                                CortexVertex cvSink = TraversalUtils.findVertex(g, sink);

                                if (cvSource != null) {
                                    if (cvSink == null) {
                                        List<CortexVertex> cvs = TraversalUtils.toWalk(g, source, c);
                                        contig = TraversalUtils.toContig(cvs);

                                        int largestDovetail = 0;

                                        for (int dovetail = 1; dovetail < GRAPH.getKmerSize(); dovetail++) {
                                            String tail = contig.substring(contig.length() - dovetail, contig.length());
                                            String head = sink.substring(0, dovetail);

                                            if (tail.equals(head) && dovetail > largestDovetail) {
                                                largestDovetail = dovetail;
                                                finalOverlap = head;
                                            }
                                        }

                                        if (largestDovetail >= 10 && largestDovetail > 0) {
                                            contig += sink.substring(largestDovetail, sink.length());
                                            finalDovetail = largestDovetail;
                                            finalColor = c;
                                            offset = left;
                                        }
                                    } else {
                                        PathFinder pf = new PathFinder(g, c);
                                        List<GraphPath<CortexVertex, CortexEdge>> gps = pf.getPaths(cvSource, cvSink);

                                        for (GraphPath<CortexVertex, CortexEdge> gp : gps) {
                                            List<CortexVertex> cvs = gp.getVertexList();
                                            contig = TraversalUtils.toContig(cvs);
                                            finalDovetail = contig.length();
                                            finalOverlap = contig;
                                            finalColor = c;
                                            offset = left;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    log.info("  region: {}-{}", region.getFirst(), region.getSecond());
                    if (finalColor >= 0) {
                        log.info("  - {} {} {} {} {} {}",
                                String.format("%3d", finalColor),
                                String.format("%-20s", GRAPH.getSampleName(finalColor)),
                                finalDovetail,
                                finalOverlap.length(),
                                contig.length(),
                                contig
                        );

                        SmithWaterman sw = new SmithWaterman();
                        String[] al = sw.getAlignment(childContig, contig);

                        log.info("  - kid {}", al[0]);
                        log.info("  - par {}", al[1]);

                        List<Triple<Integer, String, String>> calls = call(al);

                        for (Triple<Integer, String, String> call : calls) {
                            vout.println(Joiner.on("\t").join(contigIndex, call.getLeft() + offset, call.getMiddle(), call.getRight()));
                        }

                        log.info("");
                    }
                }
            }
        }
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

    private void printContigMatrix(int contigIndex, Set<CanonicalKmer> rois, List<CortexVertex> ws, double[][][] a, List<String> stateList, List<String> stateSeq, String suffix) {
        Map<Integer, String> stateMap = new HashMap<>();
        for (int i = 0; i < stateList.size(); i++) {
            stateMap.put(i, stateList.get(i));
        }

        try {
            PrintStream fout = new PrintStream(out.getAbsolutePath() + "/contig" + contigIndex + "." + suffix + ".matrix.txt");

            for (int i = 0; i < a.length; i++) {
                for (int s0 = 0; s0 < a[i].length; s0++) {
                    for (int s1 = 0; s1 < a[i][s0].length; s1++) {
                        String n0 = stateMap.get(s0);
                        String n1 = stateMap.get(s1);

                        fout.println(Joiner.on("\t").join(i, s0, s1, n0, n1, a[i][s0][s1], rois.contains(ws.get(i).getCanonicalKmer())));
                    }
                }
            }

            fout.close();

            PrintStream pout = new PrintStream(out.getAbsolutePath() + "/contig" + contigIndex + "." + suffix + ".path.txt");

            int sampleColor = GRAPH.getColorForSampleName(ROIS.getSampleName(0));
            for (int i = 0; i < a.length; i++) {
                int covChild = ws.get(i).getCortexRecord().getCoverage(sampleColor);

                int covMother = 0;
                for (int c : GRAPH.getColorsForSampleNames(MOTHER)) {
                    covMother = Math.max(covMother, ws.get(i).getCortexRecord().getCoverage(c));
                }

                int covFather = 0;
                for (int c : GRAPH.getColorsForSampleNames(FATHER)) {
                    covFather = Math.max(covFather, ws.get(i).getCortexRecord().getCoverage(c));
                }

                pout.println(Joiner.on("\t").join(i, stateList.indexOf(stateSeq.get(i)), stateSeq.get(i), ws.get(i).getKmerAsString(), covChild, covMother, covFather));
            }

            pout.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
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
