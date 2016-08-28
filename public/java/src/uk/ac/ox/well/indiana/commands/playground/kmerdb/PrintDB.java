package uk.ac.ox.well.indiana.commands.playground.kmerdb;

import com.google.common.base.Joiner;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.HTreeMap;
import org.mapdb.Serializer;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.File;
import java.util.Map;
import java.util.TreeMap;

public class PrintDB extends Module {
    @Argument(fullName="db", shortName="db", doc="Novel kmer db")
    public File DB_FILE;

    @Override
    public void execute() {
        DB db = DBMaker
                .fileDB(DB_FILE)
                .readOnly()
                .fileMmapEnable()
                .make();

        int kmerSize = db.atomicInteger("kmerSize").open().get();
        int kmerBits = db.atomicInteger("kmerBits").open().get();

        HTreeMap<long[], Map<String, Object>> nkdb = db.hashMap("novelKmers")
                .keySerializer(Serializer.LONG_ARRAY)
                .valueSerializer(Serializer.JAVA)
                .open();

        long numNovelRecords = 0;
        for (long[] bk : nkdb.getKeys()) {
            String sk = new String(CortexUtils.decodeBinaryKmer(bk, kmerSize, kmerBits));

            Map<String, Object> m = nkdb.get(bk);

            log.info("{}, '{}'", sk, Joiner.on(";").withKeyValueSeparator("=").join(m));

            numNovelRecords++;
        }

        log.info("Number of novel records: {}, {} {}", numNovelRecords, kmerSize, kmerBits);
    }
}
