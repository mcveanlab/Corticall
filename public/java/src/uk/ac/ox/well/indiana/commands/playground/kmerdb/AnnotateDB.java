package uk.ac.ox.well.indiana.commands.playground.kmerdb;

import com.google.common.base.Joiner;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.HTreeMap;
import org.mapdb.Serializer;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

public class AnnotateDB extends Module {
    @Argument(fullName="db", shortName="db", doc="Novel kmer db")
    public File DB_FILE;

    @Argument(fullName="dirty", shortName="d", doc="Dirty graph", required=false)
    public CortexGraph DIRTY;

    @Argument(fullName="clean", shortName="c", doc="Clean graph", required=false)
    public CortexGraph CLEAN;

    @Override
    public void execute() {
        DB db = DBMaker
                .fileDB(DB_FILE)
                .fileMmapEnable()
                .make();

        int kmerSize = db.atomicInteger("kmerSize").open().get();
        int kmerBits = db.atomicInteger("kmerBits").open().get();

        HTreeMap<long[], HashMap<String, Object>> nkdb = db.hashMap("novelKmers")
                .keySerializer(Serializer.LONG_ARRAY)
                .valueSerializer(Serializer.JAVA)
                .counterEnable()
                .open();

        ProgressMeter pm = new ProgressMeterFactory()
                .maxRecord(nkdb.size())
                .header("Processing database...")
                .message("records processed")
                .make(log);

        for (long[] bk : nkdb.getKeys()) {
            String sk = new String(CortexUtils.decodeBinaryKmer(bk, kmerSize, kmerBits));

            HashMap<String, Object> m = nkdb.get(bk);

            m.put("coverage_dirty", "0");

            if (DIRTY != null) {
                CortexRecord dr = DIRTY.findRecord(sk);
                if (dr != null) {
                    int cov = dr.getCoverage(0);

                    m.put("coverage_dirty", String.valueOf(cov));
                }
            }

            m.put("coverage_clean", "0");

            if (CLEAN != null) {
                CortexRecord cr = CLEAN.findRecord(sk);
                if (cr != null) {
                    int cov = cr.getCoverage(0);

                    m.put("coverage_clean", String.valueOf(cov));
                }
            }

            nkdb.put(bk, m);

            pm.update();
        }

        db.commit();
        nkdb.close();
        db.close();
    }
}
