package uk.ac.ox.well.indiana.commands.playground.index.links;

import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.HTreeMap;
import org.mapdb.Serializer;
import org.mapdb.serializer.SerializerArrayTuple;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;

import java.io.File;
import java.util.ArrayList;
import java.util.Map;
import java.util.NavigableSet;
import java.util.TreeMap;

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
        //storeLinks(db, sourceMap);

        for (CortexRecord cr : GRAPH) {
            for (int c = 0; c < LINKS.size(); c++) {
                CortexLinks cl = LINKS.get(c);
                CortexLinksRecord clr = cl.get(cr.getCortexKmer());

                for (CortexJunctionsRecord cjr : clr.getJunctions()) {
                    log.info("{} {}", c, cjr);
                }
            }
        }

        db.close();
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
