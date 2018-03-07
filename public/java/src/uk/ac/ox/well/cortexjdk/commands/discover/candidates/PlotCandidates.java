package uk.ac.ox.well.cortexjdk.commands.discover.candidates;

import com.google.common.base.Joiner;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.HTreeMap;
import org.mapdb.Serializer;
import org.mapdb.serializer.SerializerArray;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.io.*;
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

    @Argument(fullName="mother", shortName="m", doc="Mother's sample name")
    public LinkedHashSet<String> MOTHER;

    @Argument(fullName="father", shortName="f", doc="Father's sample name")
    public LinkedHashSet<String> FATHER;

    @Output
    public PrintStream out;

    @Output(fullName="rout1", shortName="ro1", doc="Reads out end 1")
    public PrintStream rout1;

    @Output(fullName="rout2", shortName="ro2", doc="Reads out end 2")
    public PrintStream rout2;

    @Output(fullName="bout", shortName="bo", doc="Bed")
    public PrintStream bout;

    @Output(fullName="fmout", shortName="fmo", doc="Fasta mother out")
    public PrintStream fmout;

    @Output(fullName="ffout", shortName="ffo", doc="Fasta father out")
    public PrintStream ffout;

    @Output(fullName="rmout1", shortName="rmo1", doc="Reads out end 1")
    public PrintStream rmout1;

    @Output(fullName="rmout2", shortName="rmo2", doc="Reads out end 2")
    public PrintStream rmout2;

    @Output(fullName="rfout1", shortName="rfo1", doc="Reads out end 1")
    public PrintStream rfout1;

    @Output(fullName="rfout2", shortName="rfo2", doc="Reads out end 2")
    public PrintStream rfout2;

    @Override
    public void execute() {
        log.info("Loading contigs from database...");
        DB dbIn = DBMaker.fileDB(DB_FILE).readOnly().fileMmapEnable().make();

        HTreeMap<Integer, CortexVertex[]> contigsMap = initializeDbContigsMap(dbIn);
        HTreeMap<Integer, FastqRecord[]> readsMapEnd1 = initializeDbReadsMap(dbIn, true);
        HTreeMap<Integer, FastqRecord[]> readsMapEnd2 = initializeDbReadsMap(dbIn, false);

        log.info("  contigs={} readsets=[{} {}]", contigsMap.size(), readsMapEnd1.size(), readsMapEnd2.size());

        int sampleColor = GRAPH.getColorForSampleName(SAMPLE);

        Set<CanonicalKmer> rois = loadRois(ROIS);

        Set<Integer> contigIndices = new TreeSet<>(contigsMap.getKeys());
        for (Integer contigIndex : contigIndices) {
            List<CortexVertex> w = getWalk(contigsMap, contigIndex);
            String childContig = TraversalUtils.toContig(w);

            Map<CanonicalKmer, Integer> positions = new HashMap<>();
            for (int i = 0; i < w.size(); i++) {
                positions.put(w.get(i).getCanonicalKmer(), i);
            }

            out.println(">contig" + contigIndex);
            out.println(childContig);

            List<FastqRecord> readsEnd1 = getReads(readsMapEnd1, contigIndex);
            List<FastqRecord> readsEnd2 = getReads(readsMapEnd2, contigIndex);

            for (int i = 0; i < readsEnd1.size(); i++) {
                rout1.println(readsEnd1.get(i));
                rout2.println(readsEnd2.get(i));
            }

            List<Pair<Integer, Integer>> regions = getRegions(rois, w);
            IntervalTreeMap<Interval> itm = new IntervalTreeMap<>();
            for (Pair<Integer, Integer> region : regions) {
                Interval it = new Interval("contig" + contigIndex, region.getFirst(), region.getSecond());
                itm.put(it, it);
            }

            log.info("{} {} {}", contigIndex, childContig.length(), TraversalUtils.toContig(w));
            for (Pair<Integer, Integer> region : regions) {
                log.info("  region: {} {}", region.getFirst(), region.getSecond());
                bout.println(Joiner.on("\t").join("contig" + contigIndex, region.getFirst(), region.getSecond(), "+"));
            }

            int subcontigStart = regions.get(0).getFirst() - 200;
            if (subcontigStart < 0) { subcontigStart = 0; }

            int subcontigStop  = regions.get(regions.size() - 1).getSecond() + 200;
            if (subcontigStop >= w.size()) { subcontigStop = w.size() - 1; }

            for (Set<String> background : Arrays.asList(MOTHER, FATHER)) {
                log.info("  background: {}", Joiner.on(", ").join(background));

                TraversalEngine ef = new TraversalEngineFactory()
                        .traversalColors(GRAPH.getColorsForSampleNames(background))
                        .traversalDirection(FORWARD)
                        .combinationOperator(OR)
                        .stoppingRule(ContigStopper.class)
                        .graph(GRAPH)
                        .links(LINKS)
                        .make();

                TraversalEngine er = new TraversalEngineFactory()
                        .traversalColors(GRAPH.getColorsForSampleNames(background))
                        .traversalDirection(REVERSE)
                        .combinationOperator(OR)
                        .stoppingRule(ContigStopper.class)
                        .graph(GRAPH)
                        .links(LINKS)
                        .make();

                //DirectedWeightedPseudograph<CortexVertex, CortexEdge> d = new DirectedWeightedPseudograph<>(CortexEdge.class);
                Set<CortexVertex> seen = new HashSet<>();

                Map<String, List<CortexVertex>> pContigs = new LinkedHashMap<>();
                for (int i = subcontigStart, j = subcontigStop; i <= j; i++, j--) {
                    for (int c : GRAPH.getColorsForSampleNames(background)) {
                        if (w.get(i).getCortexRecord().getCoverage(c) > 0 && !seen.contains(w.get(i))) {
                            List<CortexVertex> pw = ef.walk(w.get(i).getKmerAsString());
                            String pc = TraversalUtils.toContig(pw);

                            if (pc.length() > 0) {
                                pContigs.put(pc, pw);
                            }

                            seen.addAll(pw);
                        }

                        if (w.get(j).getCortexRecord().getCoverage(c) > 0 && !seen.contains(w.get(j))) {
                            List<CortexVertex> pw = er.walk(w.get(j).getKmerAsString());
                            String pc = TraversalUtils.toContig(pw);

                            if (pc.length() > 0) {
                                pContigs.put(pc, pw);
                            }

                            seen.addAll(pw);
                        }
                    }
                }

                List<String> contigList = new ArrayList<>(pContigs.keySet());
                for (int i = 0; i < contigList.size(); i++) {
                    if (background.equals(MOTHER)) {
                        fmout.println(">contig" + contigIndex + "_s" + i);
                        fmout.println(contigList.get(i));
                        fmout.flush();
                    } else if (background.equals(FATHER)) {
                        ffout.println(">contig" + contigIndex + "_s" + i);
                        ffout.println(contigList.get(i));
                        ffout.flush();
                    }
                }

                Map<String, Map<CanonicalKmer, Integer>> contigPositions = new HashMap<>();
                for (String pContig : pContigs.keySet()) {
                    Map<CanonicalKmer, Integer> cp = new HashMap<>();
                    for (int i = 0; i < pContigs.get(pContig).size(); i++) {
                        cp.put(pContigs.get(pContig).get(i).getCanonicalKmer(), i);
                    }

                    contigPositions.put(pContig, cp);
                }

                for (int i = 0; i < readsEnd1.size(); i++) {
                    FastqRecord fend1 = readsEnd1.get(i);
                    FastqRecord fend2 = readsEnd2.get(i);

                    Pair<Integer, Integer> cExtent = getExtent(positions, fend1, fend2);

                    Interval it = new Interval("contig" + contigIndex, cExtent.getFirst(), cExtent.getSecond());

                    Collection<Interval> ol = itm.getOverlapping(it);

                    if (ol.size() > 0) {
                        if (background.equals(MOTHER)) {
                            rmout1.println(fend1);
                            rmout2.println(fend2);
                            rmout1.flush();
                            rmout2.flush();
                        } else if (background.equals(FATHER)) {
                            rfout1.println(fend1);
                            rfout2.println(fend2);
                            rfout1.flush();
                            rfout2.flush();
                        }

                        Set<Interval> pExtents = new TreeSet<>(Comparator.comparingInt(Interval::length).reversed());
                        for (String pContig : contigPositions.keySet()) {
                            Pair<Integer, Integer> pExtent1 = getExtent(contigPositions.get(pContig), fend1);
                            Pair<Integer, Integer> pExtent2 = getExtent(contigPositions.get(pContig), fend2);
                            Pair<Integer, Integer> pExtent = getExtent(contigPositions.get(pContig), fend1, fend2);

                            if (pExtent.getSecond() > pExtent.getFirst() && pExtent.getFirst() >= 0 && !pExtent.equals(pExtent1) && !pExtent.equals(pExtent2)) {
                                Interval pit = new Interval("contig" + contigIndex, pExtent.getFirst(), pExtent.getSecond());
                                pExtents.add(pit);
                            }
                        }

                        for (Interval pit : pExtents) {
                            log.info("  read {} ({}): {}:{}-{} ({}) {}:{}-{} ({}) {} {} {}", i, background.iterator().next(), it.getContig(), it.getStart(), it.getEnd(), it.length(), pit.getContig(), pit.getStart(), pit.getEnd(), pit.length(), fend1.getReadHeader(), fend1.getReadString(), fend2.getReadString());
                        }

                        //log.info("  {} {}", it.length(), pExtents.size() == 0 ? 0 : pExtents.iterator().next().length());
                    }
                }
            }
        }
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
