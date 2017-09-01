package uk.ac.ox.well.cortexjdk.commands.index.links;

import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.HTreeMap;
import org.mapdb.Serializer;
import org.mapdb.serializer.SerializerArrayTuple;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinksIterable;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;

import java.io.File;
import java.util.NavigableSet;

public class IndexLinks extends Module {
    @Argument(fullName="links", shortName="l", doc="Links")
    public CortexLinksIterable LINKS;

    @Argument(fullName="source", shortName="s", doc="Link source")
    public String SOURCE;

    @Override
    public void execute() {
        File dbFile = new File(LINKS.getCortexLinksFile().getAbsolutePath() + ".k" + LINKS.getKmerSize() + ".linkdb");

        DB db = DBMaker
                .fileDB(dbFile)
                .fileMmapEnable()
                .make();

        storeVersion(db);
        storeSourceTable(db);
        storeColorTable(db);
        storeLinks(db);

        db.close();
    }

    private void storeVersion(DB db) {
        db.atomicInteger("version", 1).create();

        db.commit();
    }

    private void storeSourceTable(DB db) {
        HTreeMap<Integer, String> sources = db.hashMap("sources")
                .keySerializer(Serializer.INTEGER)
                .valueSerializer(Serializer.STRING)
                .create();

        sources.put(0, SOURCE);

        db.commit();
    }

    private void storeColorTable(DB db) {
        HTreeMap<Integer, String> header = db.hashMap("colors")
                .keySerializer(Serializer.INTEGER)
                .valueSerializer(Serializer.STRING)
                .create();

        for (int c = 0; c < LINKS.getNumColors(); c++) {
            String sampleName = LINKS.getColor(c).getSampleName();
            header.put(c, sampleName);
        }

        db.commit();
    }

    private void storeLinks(DB db) {
        NavigableSet<Object[]> linkIndex = db.treeSet("links")
                .serializer(new SerializerArrayTuple(Serializer.BYTE_ARRAY, Serializer.BYTE_ARRAY, Serializer.BYTE_ARRAY))
                .create();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing links")
                .message("links processed")
                .maxRecord(LINKS.getNumKmersWithLinks())
                .make(log);

        int numRecords = 0, numLinks = 0;
        for (CortexLinksRecord clr : LINKS) {
            for (CortexJunctionsRecord cjr : clr.getJunctions()) {
                linkIndex.add(new Object[]{clr.getKmer().getKmerAsBytes(), cjr.toString().getBytes(), "0".getBytes()});
                numLinks++;
            }
            numRecords++;

            pm.update();
        }

        log.info("Wrote {} links in {} records", numLinks, numRecords);

        db.commit();
    }
}
