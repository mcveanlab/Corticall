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
import uk.ac.ox.well.cortexjdk.utils.alignment.mosaic.MosaicAligner;
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

    @Output(fullName="targetOut", shortName="to", doc="Target out")
    public PrintStream tout;

    @Output(fullName="panelOut", shortName="po", doc="Panel out")
    public PrintStream pout;

    @Override
    public void execute() {
        log.info("Loading contigs from database...");
        DB dbIn = DBMaker.fileDB(DB_FILE).readOnly().fileMmapEnable().make();

        HTreeMap<Integer, CortexVertex[]> contigsMap = initializeDbContigsMap(dbIn);
        //HTreeMap<Integer, FastqRecord[]> readsMapEnd1 = initializeDbReadsMap(dbIn, true);
        //HTreeMap<Integer, FastqRecord[]> readsMapEnd2 = initializeDbReadsMap(dbIn, false);

        String sample = ROIS.getSampleName(0);
        Set<CanonicalKmer> rois = loadRois(ROIS);

        MosaicAligner ma = new MosaicAligner();

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

            Map<String, String> targets = new LinkedHashMap<>();
            List<Pair<String, Interval>> ls = bf.getStates(PARENT, 10);
            for (Pair<String, Interval> l : ls) {
                Interval it = l.getSecond();
                if (it.getStart() < 0) {
                    it = new Interval(it.getContig(), 1, it.getEnd(), it.isNegativeStrand(), it.getName());
                }

                int maxLength = BACKGROUNDS.get(l.getFirst()).getReferenceSequence().getSequenceDictionary().getSequence(it.getContig()).getSequenceLength();
                if (it.getEnd() >= maxLength) {
                    it = new Interval(it.getContig(), it.getStart(), maxLength - 1, it.isNegativeStrand(), it.getName());
                }

                String seq = BACKGROUNDS.get(l.getFirst()).find(it);
                String targetName = l.getFirst() + ":" + it.getContig() + ":" + it.getStart() + "-" + it.getEnd() + ":" + (it.isPositiveStrand() ? "+" : "-");
                targets.put(targetName, seq);
            }

            ma.align(TraversalUtils.toContig(ws), targets);

            log.info("\n{}", ma);

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
