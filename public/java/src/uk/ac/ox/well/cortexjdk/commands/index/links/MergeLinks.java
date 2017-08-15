package uk.ac.ox.well.cortexjdk.commands.index.links;

import com.google.common.base.Joiner;
import org.apache.commons.math3.util.Pair;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.HTreeMap;
import org.mapdb.Serializer;
import org.mapdb.serializer.SerializerArrayTuple;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;

import java.io.File;
import java.util.*;

public class MergeLinks extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Links")
    public ArrayList<CortexLinks> LINKS;

    @Output
    public File out;

    @Override
    public void execute() {
        Map<String, Integer> sourceMap = getSourceMap();

        DB db = DBMaker
                .fileDB(out)
                .fileMmapEnable()
                .make();

        storeVersion(db, 1);
        storeSourceTable(db, sourceMap);
        storeColorTable(db);
        storeLinks(db, sourceMap);

        db.close();
    }

    private void storeLinks(DB db, Map<String, Integer> sourceMap) {
        NavigableSet<Object[]> linkIndex = db.treeSet("links")
                .serializer(new SerializerArrayTuple(Serializer.BYTE_ARRAY, Serializer.BYTE_ARRAY, Serializer.BYTE_ARRAY))
                .create();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing records")
                .message("records processed")
                .maxRecord(GRAPH.getNumRecords())
                .updateRecord(GRAPH.getNumRecords() / 1000)
                .make(log);

        int numLinks = 0;
        for (CortexRecord cr : GRAPH) {
            Map<Pair<String, Boolean>, Map<Integer, Integer>> covs = new LinkedHashMap<>();
            Map<Pair<String, Boolean>, Map<Integer, Set<Integer>>> srcs = new LinkedHashMap<>();

            for (int i = 0; i < LINKS.size(); i++) {
                int color = GRAPH.getColorForSampleName(LINKS.get(i).getSampleNameForColor(0));

                CortexLinks cl = LINKS.get(i);
                CortexLinksRecord clr = cl.get(cr.getCortexKmer());

                for (CortexJunctionsRecord cjr : clr.getJunctions()) {
                    Pair<String, Boolean> key = new Pair<>(cjr.getJunctions(), cjr.isForward());

                    if (!covs.containsKey(key)) {
                        covs.put(key, new TreeMap<>());
                        srcs.put(key, new TreeMap<>());
                    }

                    if (!covs.get(key).containsKey(color)) {
                        covs.get(key).put(color, 0);
                        srcs.get(key).put(color, new TreeSet<>());
                    }

                    if (cjr.getCoverage(0) > covs.get(key).get(color)) {
                        covs.get(key).put(color, cjr.getCoverage(0));
                    }

                    srcs.get(key).get(color).addAll(getSourceIndices(sourceMap, LINKS.get(i).getSources()));
                }
            }

            for (Pair<String, Boolean> key : covs.keySet()) {
                CortexJunctionsRecord ncjr = new CortexJunctionsRecord(key.getSecond(), 0, key.getFirst().length(), asArray(covs.get(key)), key.getFirst(), asMultiArray(srcs.get(key)));

                linkIndex.add(new Object[]{cr.getKmerAsBytes(), ncjr.toString().getBytes(), ncjr.getSources().getBytes()});
                numLinks++;
            }

            pm.update("records processed, " + numLinks + " stored");
        }

        db.commit();
    }

    private int[] asArray(Map<Integer, Integer> covs) {
        int[] cova = new int[GRAPH.getNumColors()];

        for (int c = 0; c < GRAPH.getNumColors(); c++) {
            cova[c] = covs.containsKey(c) ? covs.get(c) : 0;
        }

        return cova;
    }

    private int[][] asMultiArray(Map<Integer, Set<Integer>> srcs) {
        int[][] srcsa = new int[GRAPH.getNumColors()][];

        for (int c = 0; c < GRAPH.getNumColors(); c++) {
            if (srcs.containsKey(c)) {
                Set<Integer> ss = srcs.get(c);
                int[] ssb = new int[ss.size()];

                int i = 0;
                for (int si : ss) {
                    ssb[i] = si;
                    i++;
                }

                srcsa[c] = ssb;
            } else {
                srcsa[c] = new int[0];
            }
        }

        return srcsa;
    }

    private String asString(Map<Integer, Integer> covs, String sep) {
        List<Integer> covList = new ArrayList<>();
        for (int c = 0; c < GRAPH.getNumColors(); c++) {
            covList.add(covs.containsKey(c) ? covs.get(c) : 0);
        }

        return Joiner.on(sep).join(covList);
    }

    private String asString(Map<Integer, Set<Integer>> srcs, String isep, String esep) {
        List<String> srcList = new ArrayList<>();
        for (int c = 0; c < GRAPH.getNumColors(); c++) {
            if (srcs.containsKey(c)) {
                srcList.add(Joiner.on(isep).join(srcs.get(c)));
            } else {
                srcList.add("");
            }
        }

        return Joiner.on(esep).join(srcList);
    }

    private Set<Integer> getSourceIndices(Map<String, Integer> sourceMap, Collection<String> sources) {
        Set<Integer> sourceIndices = new TreeSet<>();

        for (String source : sources) {
            sourceIndices.add(sourceMap.get(source));
        }

        return sourceIndices;
    }

    private Map<String, Integer> getSourceMap() {
        Map<String, Integer> sourceMap = new TreeMap<>();
        int sourceId = 0;
        for (CortexLinks l : LINKS) {
            for (String source : l.getSources()) {
                if (!sourceMap.containsKey(source)) {
                    sourceMap.put(source, sourceId);
                    sourceId++;
                }
            }
        }

        return sourceMap;
    }

    private void storeVersion(DB db, int version) {
        db.atomicInteger("version", version).create();

        db.commit();
    }

    private void storeSourceTable(DB db, Map<String, Integer> sourceMap) {
        HTreeMap<Integer, String> sources = db.hashMap("sources")
                .keySerializer(Serializer.INTEGER)
                .valueSerializer(Serializer.STRING)
                .create();

        for (String source : sourceMap.keySet()) {
            sources.put(sourceMap.get(source), source);
        }

        db.commit();
    }

    private void storeColorTable(DB db) {
        HTreeMap<Integer, String> header = db.hashMap("colors")
                .keySerializer(Serializer.INTEGER)
                .valueSerializer(Serializer.STRING)
                .create();

        for (int c = 0; c < GRAPH.getNumColors(); c++) {
            String sampleName = GRAPH.getColor(c).getSampleName();
            header.put(c, sampleName);
        }

        db.commit();
    }

    /*
    private void storeLinks(DB db, Map<String, Integer> sourceMap) {
        NavigableSet<Object[]> linkIndex = db.treeSet("links")
                .serializer(new SerializerArrayTuple(Serializer.BYTE_ARRAY, Serializer.BYTE_ARRAY, Serializer.INTEGER))
                .create();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Merging links")
                .message("files processed")
                .maxRecord(LINKS.size())
                .make(log);

        for (int c = 0; c < LINKS.size(); c++) {
            for (CortexKmer ck : LINKS.get(c).keySet()) {
                for (CortexJunctionsRecord cjr : LINKS.get(c).get(ck).getJunctions()) {
                    linkIndex.add(new Object[]{ck.getKmerAsBytes(), cjr.toString().getBytes(), LINKS.get(c).getSourceForIndex(0)});
                }
            }

            pm.update();
        }

        //log.info("Wrote {} links in {} records", numLinks, numRecords);

        db.commit();
    }
    */
}
