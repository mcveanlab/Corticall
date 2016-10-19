package uk.ac.ox.well.indiana.commands.playground.kmerdb;

import com.google.common.base.Joiner;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.ReferenceSequence;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.HTreeMap;
import org.mapdb.Serializer;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.commands.playground.index.KmerIndex;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.BlastAligner;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.DustMasker;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.util.*;

public class AnnotateDbWithBlast extends Module {
    @Argument(fullName="db", shortName="db", doc="Novel kmer db")
    public File DB_FILE;

    @Output
    public File out;

    @Override
    public void execute() {
        DB dbi = DBMaker
                 .fileDB(DB_FILE)
                 .fileMmapEnable()
                 .readOnly()
                 .make();

        int kmerSize = dbi.atomicInteger("kmerSize").open().get();
        int kmerBits = dbi.atomicInteger("kmerBits").open().get();

        HTreeMap<long[], HashMap<String, Object>> nkdb = dbi.hashMap("novelKmers")
                .keySerializer(Serializer.LONG_ARRAY)
                .valueSerializer(Serializer.JAVA)
                .counterEnable()
                .open();

        ProgressMeter pm = new ProgressMeterFactory()
                .maxRecord(nkdb.size())
                .updateRecord(nkdb.size()/100)
                .header("Processing database...")
                .message("records processed")
                .make(log);

        Map<Integer, ReferenceSequence> contigs = new HashMap<>();

        for (long[] bk : nkdb.getKeys()) {
            HashMap<String, Object> m = nkdb.get(bk);

            String contig = (String) m.get("contig");
            int group = (Integer) m.get("group");

            ReferenceSequence rseq = new ReferenceSequence(String.valueOf(group), group, contig.getBytes());
            contigs.put(group, rseq);

            pm.update();
        }

        BlastAligner ba = new BlastAligner();
        Map<Integer, Map<String, Object>> ls = ba.align(contigs.values());

        DustMasker dm = new DustMasker();
        Map<Integer, Map<String, Object>> ms = dm.mask(contigs.values());

        if (out.exists()) { out.delete(); }

        DB dbo = DBMaker
                 .fileDB(out)
                 .fileMmapEnable()
                 .make();

        HTreeMap<long[], HashMap<String, Object>> anndb = dbo.hashMap("novelKmers")
                .keySerializer(Serializer.LONG_ARRAY)
                .valueSerializer(Serializer.JAVA)
                .create();

        dbo.atomicInteger("kmerSize", kmerSize).create();
        dbo.atomicInteger("kmerBits", kmerBits).create();

        pm.reset();

        for (long[] bk : nkdb.getKeys()) {
            HashMap<String, Object> m = nkdb.get(bk);

            int group = (Integer) m.get("group");

            if (ls.containsKey(group)) {
                m.putAll(ls.get(group));
            }

            m.putAll(ms.get(group));

            anndb.put(bk, m);

            pm.update();
        }

        dbo.commit();
        dbo.close();
        dbi.close();
        nkdb.close();
    }
}
