package uk.ac.ox.well.cortexjdk.commands.call.call;

import com.google.common.base.Joiner;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.HTreeMap;
import org.mapdb.Serializer;
import org.mapdb.serializer.SerializerArray;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.graph.PairedEndAlignmentInfo;
import uk.ac.ox.well.cortexjdk.utils.alignment.graph.ReadToGraphAligner;
import uk.ac.ox.well.cortexjdk.utils.alignment.graph.SingleEndAlignmentInfo;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.DestinationStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ExplorationStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.REVERSE;

public class PlotCandidates extends Module {
    @Argument(fullName = "graph", shortName = "g", doc = "Graph")
    public CortexGraph GRAPH;

    @Argument(fullName = "links", shortName = "l", doc = "Links", required=false)
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName = "rois", shortName = "r", doc = "Rois")
    public CortexGraph ROIS;

    @Argument(fullName="db", shortName="d", doc="Database")
    public File DB_FILE;

    @Argument(fullName="sample", shortName="s", doc="Sample name")
    public String SAMPLE;

    @Argument(fullName="background", shortName="b", doc="Background name")
    public ArrayList<String> BACKGROUNDS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        log.info("Loading contigs from database...");
        DB dbIn = DBMaker.fileDB(DB_FILE).readOnly().fileMmapEnable().make();

        HTreeMap<Integer, CortexVertex[]> contigsMap = initializeDbContigsMap(dbIn);
        HTreeMap<Integer, FastqRecord[]> readsMapEnd1 = initializeDbReadsMap(dbIn, true);
        HTreeMap<Integer, FastqRecord[]> readsMapEnd2 = initializeDbReadsMap(dbIn, false);

        log.info("  contigs={} readsets=[{} {}]", contigsMap.size(), readsMapEnd1.size(), readsMapEnd2.size());

        int sampleColor = GRAPH.getColorForSampleName(SAMPLE);
        Set<Integer> backgroundColors = new TreeSet<>(GRAPH.getColorsForSampleNames(BACKGROUNDS));
        Set<Integer> colors = new TreeSet<>();
        colors.add(sampleColor);
        colors.addAll(backgroundColors);
        Map<Integer, Integer> recruitmentColors = new HashMap<>();
        recruitmentColors.put(0, 22);
        recruitmentColors.put(1, 23);

        Set<CanonicalKmer> rois = loadRois(ROIS);

        for (Integer contigIndex : contigsMap.getKeys()) {
            List<CortexVertex> w = getWalk(contigsMap, contigIndex);

            List<Pair<Integer, Integer>> regions = getRegions(rois, w, 200);
            int subcontigStart = regions.get(0).getFirst();
            int subcontigStop  = regions.get(regions.size() - 1).getSecond();

            List<CortexVertex> ws = new ArrayList<>();
            for (int i = subcontigStart; i < subcontigStop; i++) {
                ws.add(w.get(i));
            }

            DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = new DirectedWeightedPseudograph<>(CortexEdge.class);
            for (int i = 0; i < ws.size(); i++) {
                g.addVertex(ws.get(i));
            }

            log.info("{} {} {} {}", contigIndex, w.size(), subcontigStart, subcontigStop);

            for (int c : colors) {
                for (int i = 0; i < ws.size() - 1; i++) {
                    CortexVertex v0 = ws.get(i);
                    CortexVertex v1 = ws.get(1);

                    Set<String> nks = convertToStrings(TraversalUtils.getAllNextKmers(v0.getCortexRecord(), !v0.getKmerAsString().equals(v0.getCanonicalKmer().getKmerAsString())).get(c));

                    if (nks.contains(v1.getKmerAsString())) {
                        g.addEdge(v0, v1, new CortexEdge(v0, v1, c, 1.0));
                    }
                }

                TraversalEngine ef = new TraversalEngineFactory()
                        .traversalColor(c)
                        .traversalDirection(FORWARD)
                        .combinationOperator(OR)
                        .stoppingRule(ExplorationStopper.class)
                        .maxBranchLength(1000)
                        .graph(GRAPH)
                        .links(LINKS)
                        .make();

                TraversalEngine er = new TraversalEngineFactory()
                        .traversalColor(c)
                        .traversalDirection(REVERSE)
                        .combinationOperator(OR)
                        .stoppingRule(ExplorationStopper.class)
                        .maxBranchLength(1000)
                        .graph(GRAPH)
                        .links(LINKS)
                        .make();

                for (int i = subcontigStart; i < subcontigStop - 1; i++) {
                    CortexVertex v0 = w.get(i);
                    CortexVertex v1 = w.get(i+1);

                    Map<Integer, Set<CortexByteKmer>> allNextEdges = TraversalUtils.getAllNextKmers(v0.getCortexRecord(), !v0.getKmerAsString().equals(v0.getCanonicalKmer().getKmerAsString()));
                    Set<String> remainingNext = convertToStrings(allNextEdges.get(sampleColor));
                    remainingNext.remove(v1.getKmerAsString());

                    if (remainingNext.size() > 0) {
                        for (String nk : remainingNext) {
                            DirectedWeightedPseudograph<CortexVertex, CortexEdge> gn = ef.dfs(nk);
                            if (gn != null) {
                                Graphs.addGraph(g, gn);
                            }
                        }
                    }

                    Map<Integer, Set<CortexByteKmer>> allPrevEdges = TraversalUtils.getAllNextKmers(v1.getCortexRecord(), !v1.getKmerAsString().equals(v1.getCanonicalKmer().getKmerAsString()));
                    Set<String> remainingPrev = convertToStrings(allPrevEdges.get(sampleColor));
                    remainingPrev.remove(v0.getKmerAsString());

                    if (remainingPrev.size() > 0) {
                        for (String pk : remainingNext) {
                            DirectedWeightedPseudograph<CortexVertex, CortexEdge> gp = er.dfs(pk);
                            if (gp != null) {
                                Graphs.addGraph(g, gp);
                            }
                        }
                    }
                }

                log.info("  -- {} {}", g.vertexSet().size(), g.edgeSet().size());
            }

            log.info("  {} {}", g.vertexSet().size(), g.edgeSet().size());

            g.

            //List<FastqRecord> readsEnd1 = getReads(readsMapEnd1, contigIndex);
            //List<FastqRecord> readsEnd2 = getReads(readsMapEnd2, contigIndex);

            /*
            Map<CanonicalKmer, Integer> contigKmers = getRegionKmers(cvs, new Pair<>(0, cvs.size() - 1));

            List<Pair<Integer, Integer>> regions = getRegions(rois, cvs);

            Set<String> lines = new TreeSet<>();

            for (int ri = 0; ri < regions.size(); ri++) {
                Pair<Integer, Integer> region = regions.get(ri);

                for (int limit : Arrays.asList(region.getFirst(), region.getSecond())) {
                    List<String> fields = new ArrayList<>();
                    fields.add("limit");
                    fields.add(String.valueOf(contigIndex));
                    fields.add(String.valueOf(limit));
                    fields.add("NA");
                    for (int bc : backgroundColors) {
                        fields.add("NA");
                        fields.add("NA");
                    }

                    lines.add(Joiner.on("\t").join(fields));
                }

                Map<CanonicalKmer, Integer> regionKmers = getRegionKmers(cvs, region);
                List<Pair<FastqRecord, FastqRecord>> regionReads = getRegionReads(readsEnd1, readsEnd2, regionKmers);

                ProgressMeter pm = new ProgressMeterFactory()
                        .header("Processing reads for contig " + contigIndex + " (len=" + cvs.size() + ") region " + ri + " (" + region.getFirst() + "-" + region.getSecond() + ")")
                        .message("reads processed")
                        .maxRecord(regionReads.size())
                        .make(log);

                ReadToGraphAligner ga = new ReadToGraphAligner(GRAPH, LINKS, colors);

                for (CanonicalKmer ck : regionKmers.keySet()) {
                    if (rois.contains(ck)) {
                        List<String> fields = new ArrayList<>();
                        fields.add("novel");
                        fields.add(String.valueOf(contigIndex));
                        fields.add(String.valueOf(regionKmers.get(ck)));
                        fields.add("NA");
                        for (int bc : backgroundColors) {
                            fields.add("NA");
                            fields.add("NA");
                        }

                        lines.add(Joiner.on("\t").join(fields));
                    }
                }

                for (int i = 0; i < regionReads.size(); i++) {
                    pm.update();

                    FastqRecord fqEnd1 = regionReads.get(i).getFirst();
                    FastqRecord fqEnd2 = regionReads.get(i).getSecond();

                    PairedEndAlignmentInfo pai = ga.align(fqEnd1, fqEnd2);

                    Pair<Integer, Integer> ext1 = getExtent(contigKmers, fqEnd1);
                    Pair<Integer, Integer> ext2 = getExtent(contigKmers, fqEnd2);
                    Pair<Integer, Integer> extent = getExtent(contigKmers, fqEnd1, fqEnd2);

                    //log.info("  {}/{} {} {} {} {} {}", i, regionReads.size(), extent.getFirst(), extent.getSecond(), Joiner.on(" ").withKeyValueSeparator(" => ").join(pai.getAlignmentDistanceMap()), fqEnd1.getReadString(), fqEnd2.getReadString());
                    //log.info("  * {} {} {}", ext1, ext2, extent);

                    List<String> fields = new ArrayList<>();
                    fields.add("read");
                    fields.add(String.valueOf(contigIndex));
                    fields.add(String.valueOf(extent.getFirst()));
                    fields.add(pai.getAlignmentDistanceMap().get(sampleColor).size() == 0 ? "NA" : String.valueOf(pai.getAlignmentDistanceMap().get(sampleColor).iterator().next()));

                    for (int bc : backgroundColors) {
                        fields.add(pai.getAlignmentDistanceMap().get(bc).size() == 0 ? "NA" : String.valueOf(pai.getAlignmentDistanceMap().get(bc).iterator().next()));
                    }

                    for (SingleEndAlignmentInfo sai : Arrays.asList(pai.getFirstEndAlignment(), pai.getSecondEndAlignment())) {
                        for (int c : backgroundColors) {
                            fields.add(String.valueOf(sai.isSplitRead(sampleColor, c)));
                        }
                    }

                    lines.add(Joiner.on("\t").join(fields));
                }
            }

            for (String line : lines) {
                out.println(line);
            }

            out.flush();
            */
        }
    }

    private Set<String> convertToStrings(Set<CortexByteKmer> bs) {
        Set<String> es = new HashSet<>();

        for (CortexByteKmer b : bs) {
            es.add(new String(b.getKmer()));
        }

        return es;
    }

    private Pair<Integer, Integer> getExtent(Map<CanonicalKmer, Integer> contigKmers, FastqRecord... fqEnds) {
        int leftmostIndex = -1;
        int rightmostIndex = -1;

        for (FastqRecord fqr : fqEnds) {
            for (int i = 0; i <= fqr.length() - GRAPH.getKmerSize(); i++) {
                CanonicalKmer ck = new CanonicalKmer(fqr.getReadString().substring(i, i + GRAPH.getKmerSize()));

                if (contigKmers.containsKey(ck)) {
                    if (leftmostIndex == -1 || contigKmers.get(ck) < leftmostIndex) {
                        leftmostIndex = contigKmers.get(ck);
                    }

                    if (rightmostIndex == -1 || contigKmers.get(ck) > rightmostIndex) {
                        rightmostIndex = contigKmers.get(ck);
                    }
                }
            }
        }

        return new Pair<>(leftmostIndex, rightmostIndex);
    }

    @NotNull
    private List<Pair<Integer, Integer>> getRegions(Set<CanonicalKmer> rois, List<CortexVertex> cvs, int padding) {
        List<Pair<Integer, Integer>> regions = new ArrayList<>();
        for (int i = 0; i < cvs.size(); i++) {
            CanonicalKmer currentKmer = cvs.get(i).getCanonicalKmer();

            if (rois.contains(currentKmer)) {
                int regionStart = i - GRAPH.getKmerSize() - padding > 0 ? i - GRAPH.getKmerSize() - padding : 0;
                int regionStop = i + GRAPH.getKmerSize() + padding < cvs.size() - 1 ? i + GRAPH.getKmerSize() + padding : cvs.size() - 1;

                regions.add(new Pair<>(regionStart, regionStop));

                i = regionStop;
            }
        }
        return regions;
    }

    @NotNull
    private List<Pair<FastqRecord, FastqRecord>> getRegionReads(List<FastqRecord> readsEnd1, List<FastqRecord> readsEnd2, Map<CanonicalKmer, Integer> regionKmers) {
        List<Pair<FastqRecord, FastqRecord>> regionReads = new ArrayList<>();
        for (int i = 0; i < readsEnd1.size(); i++) {
            FastqRecord fqEnd1 = readsEnd1.get(i);
            FastqRecord fqEnd2 = readsEnd2.get(i);

            boolean addPair = false;

            for (int j = 0; j <= fqEnd1.length() - GRAPH.getKmerSize(); j++) {
                CanonicalKmer ck = new CanonicalKmer(fqEnd1.getReadString().substring(j, j + GRAPH.getKmerSize()));

                if (regionKmers.containsKey(ck)) {
                    addPair = true;
                    break;
                }
            }

            for (int j = 0; j <= fqEnd2.length() - GRAPH.getKmerSize(); j++) {
                CanonicalKmer ck = new CanonicalKmer(fqEnd2.getReadString().substring(j, j + GRAPH.getKmerSize()));

                if (regionKmers.containsKey(ck)) {
                    addPair = true;
                    break;
                }
            }

            if (addPair) {
                regionReads.add(new Pair<>(fqEnd1, fqEnd2));
            }
        }
        return regionReads;
    }

    @NotNull
    private Map<CanonicalKmer, Integer> getRegionKmers(List<CortexVertex> cvs, Pair<Integer, Integer> region) {
        Map<CanonicalKmer, Integer> regionKmers = new HashMap<>();
        for (int i = region.getFirst(); i < region.getSecond(); i++) {
            regionKmers.put(cvs.get(i).getCanonicalKmer(), i);
        }
        return regionKmers;
    }

    private List<FastqRecord> getReads(HTreeMap<Integer, FastqRecord[]> readsMap, int contigIndex) {
        List<FastqRecord> fqs = new ArrayList<>();

        Object[] oa = readsMap.get(contigIndex);

        for (Object o : oa) {
            fqs.add((FastqRecord) o);
        }

        return fqs;
    }

    private List<CortexVertex> getWalk(HTreeMap<Integer, CortexVertex[]> contigsMap, int contigIndex) {
        List<CortexVertex> cvs = new ArrayList<>();

        Object[] oa = contigsMap.get(contigIndex);

        for (Object o : oa) {
            cvs.add((CortexVertex) o);
        }

        return cvs;
    }

    private HTreeMap<Integer, FastqRecord[]> initializeDbReadsMap(DB dbIn, boolean firstEnd) {
        HTreeMap<Integer, FastqRecord[]> dbReads = dbIn
                .hashMap("readsEnd" + (firstEnd ? "1" : "2"))
                .keySerializer(Serializer.INTEGER)
                .valueSerializer(new SerializerArray(Serializer.JAVA))
                .counterEnable()
                .open();

        return dbReads;
    }

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
