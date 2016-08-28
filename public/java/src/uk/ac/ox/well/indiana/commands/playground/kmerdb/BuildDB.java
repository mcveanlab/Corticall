package uk.ac.ox.well.indiana.commands.playground.kmerdb;

import kotlin.Pair;
import org.mapdb.*;
import org.mapdb.serializer.SerializerArray;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class BuildDB extends Module {
    @Argument(fullName="dirtyGraph", shortName="d", doc="Dirty graph")
    public CortexGraph DIRTY;

    @Argument(fullName="childColor", shortName="cc", doc="Child color")
    public Integer CHILD_COLOR = 0;

    @Output
    public File out;

    @Override
    public void execute() {
        List<long[]> binaryKmers = new ArrayList<>();

        ProgressMeter pm = new ProgressMeterFactory()
                .maxRecord(DIRTY.getNumRecords())
                .header("Processing graph...")
                .message("records processed")
                .make(log);

        long numNovelRecords = 0;
        for (CortexRecord cr : DIRTY) {
            if (CortexUtils.isNovelKmer(cr, CHILD_COLOR)) {
                binaryKmers.add(cr.getBinaryKmer());

                numNovelRecords++;
            }

            pm.update("records processed (" + numNovelRecords + " novel so far)");
        }

        log.info("  {} total novel records", numNovelRecords);

        if (out.exists()) {
            out.delete();
        }

        DB db = DBMaker
                .fileDB(out)
                .fileMmapEnable()
                .make();

        Atomic.Integer aKmerSize = db.atomicInteger("kmerSize", DIRTY.getKmerSize()).create();
        Atomic.Integer aKmerBits = db.atomicInteger("kmerBits", DIRTY.getKmerBits()).create();

        HTreeMap<long[], Map<String, Object>> nkdb = db.hashMap("novelKmers")
                .keySerializer(Serializer.LONG_ARRAY)
                .valueSerializer(Serializer.JAVA)
                .create();

        pm = new ProgressMeterFactory()
                .maxRecord(binaryKmers.size())
                .header("Writing records...")
                .message("records written")
                .make(log);

        for (long[] bk : binaryKmers) {
            nkdb.put(bk, new HashMap<>());

            pm.update();
        }

        db.commit();
        nkdb.close();
        db.close();
    }
}
