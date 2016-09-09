package uk.ac.ox.well.indiana.commands.playground.kmerdb;

import com.google.common.base.Joiner;
import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.HTreeMap;
import org.mapdb.Serializer;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.commands.playground.index.KmerIndex;
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

public class AnnotateDB extends Module {
    @Argument(fullName="db", shortName="db", doc="Novel kmer db")
    public File DB_FILE;

    @Argument(fullName="dirty", shortName="d", doc="Dirty graph", required=false)
    public CortexGraph DIRTY;

    @Argument(fullName="clean", shortName="c", doc="Clean graph", required=false)
    public CortexGraph CLEAN;

    @Argument(fullName="bam", shortName="b", doc="BAM file", required=false)
    public File BAM_FILE;

    @Output
    public File out;

    @Override
    public void execute() {
        SamReader sr = null;
        KmerIndex ki = null;
        if (BAM_FILE != null) {
            sr = SamReaderFactory.make()
                    .setOption(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, true)
                    .setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, true)
                    .open(BAM_FILE);

            ki = new KmerIndex(BAM_FILE, CLEAN.getKmerSize(), false);
        }

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

        Map<long[], HashMap<String, Object>> mo = new HashMap<>();

        for (long[] bk : nkdb.getKeys()) {
            String sk = new String(CortexUtils.decodeBinaryKmer(bk, kmerSize, kmerBits));

            HashMap<String, Object> m = nkdb.get(bk);

            CortexRecord dr = null;
            if (DIRTY != null) {
                dr = DIRTY.findRecord(sk);
                if (dr != null) {
                    m.put("coverage_dirty", dr.getCoverage(0));
                    m.put("edges_dirty", dr.getEdgesAsBytes(0));
                }
            }

            CortexRecord cr = null;
            if (CLEAN != null) {
                cr = CLEAN.findRecord(sk);
                if (cr != null) {
                    m.put("coverage_dirty", cr.getCoverage(0));
                    m.put("edges_clean", cr.getEdgesAsBytes(0));
                }
            }

            if (sr != null) {
                Set<SAMRecord> srs = new HashSet<>();

                for (long[] c : ki.find(sk)) {
                    SAMFileSpan sfs = new BAMFileSpan(new Chunk(c[0], c[1]));

                    SAMRecordIterator recs = sr.indexing().iterator(sfs);

                    while (recs.hasNext()) {
                        SAMRecord s = recs.next();

                        if (s.getReadString().contains(sk) || s.getReadString().contains(SequenceUtils.reverseComplement(sk))) {
                            srs.add(s);
                        }
                    }

                    recs.close();
                }

                for (SAMRecord s : srs) {
                    String rs = s.getReadString().contains(sk) ? s.getReadString() : SequenceUtils.reverseComplement(s.getReadString());

                    log.info("{} {}\n{}\n{}\n{}{}",
                            dr == null ? 0 : dr.getCoverage(0),
                            cr == null ? 0 : cr.getCoverage(0),
                            s.getSAMString().trim(),
                            rs,
                            StringUtil.repeatCharNTimes(' ', rs.indexOf(sk)),
                            sk
                    );

                }

                m.put("reads", srs);
            }

            mo.put(bk, m);

            pm.update();
        }

        nkdb.close();
        dbi.close();

        DB dbo = DBMaker
                 .fileDB(out)
                 .fileMmapEnable()
                 .make();

        nkdb = dbo.hashMap("novelKmers")
                .keySerializer(Serializer.LONG_ARRAY)
                .valueSerializer(Serializer.JAVA)
                .create();

        for (long[] bk : mo.keySet()) {
            nkdb.put(bk, mo.get(bk));
        }

        dbo.commit();
        dbo.close();
        nkdb.close();
    }
}
