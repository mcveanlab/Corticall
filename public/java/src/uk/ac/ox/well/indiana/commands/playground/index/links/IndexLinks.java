package uk.ac.ox.well.indiana.commands.playground.index.links;

import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.Serializer;
import org.mapdb.serializer.SerializerArrayTuple;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksIterable;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;

import java.io.File;
import java.util.NavigableSet;

public class IndexLinks extends Module {
    @Argument(fullName="links", shortName="l", doc="Links")
    public CortexLinksIterable LINKS;

    @Override
    public void execute() {
        File dbFile = new File(LINKS.getCortexLinksFile().getAbsolutePath() + ".linkdb");

        DB db = DBMaker
                .fileDB(dbFile)
                .fileMmapEnable()
                .make();

        NavigableSet<Object[]> linkIndex = db.treeSet("links")
                .serializer(new SerializerArrayTuple(Serializer.BYTE_ARRAY, Serializer.BYTE_ARRAY))
                .create();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing links")
                .message("links processed")
                .maxRecord(LINKS.getNumKmersWithLinks())
                .make(log);

        int numRecords = 0, numLinks = 0;
        for (CortexLinksRecord clr : LINKS) {
            for (CortexJunctionsRecord cjr : clr.getJunctions()) {
                linkIndex.add(new Object[]{clr.getKmer().getKmerAsBytes(), cjr.toString().getBytes()});
                numLinks++;
            }
            numRecords++;

            pm.update();
        }

        db.commit();

        db.close();

        log.info("Wrote {} links in {} records", numLinks, numRecords);
    }
}
