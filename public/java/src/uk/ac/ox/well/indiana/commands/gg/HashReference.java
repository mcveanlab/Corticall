package uk.ac.ox.well.indiana.commands.gg;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;

import java.io.File;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class HashReference extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference sequence")
    public File REF_FILE;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        File dbout = new File(REF_FILE.getAbsoluteFile() + ".kmerdb");

        if (dbout.exists()) {
            throw new IndianaException("Index '" + dbout.getAbsolutePath() + "' already exists.");
        }

        if (!dbout.exists()) {
            log.info("Indexing reference to '{}'...", dbout.getAbsolutePath());

            DB db = DBMaker.newFileDB(dbout)
                    .transactionDisable()
                    .mmapFileEnable()
                    .cacheSize(1000000)
                    .closeOnJvmShutdown()
                    .make();

            Map<String, Set<CompactSerializableInterval>> index = db.getHashMap("index" + KMER_SIZE);

            FastaSequenceFile ref = new FastaSequenceFile(REF_FILE, true);
            ReferenceSequence rseq;
            while ((rseq = ref.nextSequence()) != null) {
                log.info("  {}", rseq.getName());

                String seq = new String(rseq.getBases());

                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    String sk = seq.substring(i, i + KMER_SIZE);

                    Set<CompactSerializableInterval> intervals = index.containsKey(sk) ? index.get(sk) : new HashSet<CompactSerializableInterval>();

                    intervals.add(new CompactSerializableInterval(rseq.getContigIndex(), i));

                    index.put(sk, intervals);
                }

                db.commit();
            }

            log.info("Compacting index...");
            db.compact();

            db.close();
        }

        log.info("Verifying index...");

        DB db = DBMaker.newFileDB(dbout)
                .transactionDisable()
                .mmapFileEnable()
                .cacheSize(1000000)
                .closeOnJvmShutdown()
                .readOnly()
                .make();

        Map<String, Set<CompactSerializableInterval>> newindex = db.getHashMap("index" + KMER_SIZE);

        FastaSequenceFile ref = new FastaSequenceFile(REF_FILE, true);
        ReferenceSequence rseq;
        while ((rseq = ref.nextSequence()) != null) {
            log.info("  {}", rseq.getName());

            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                String sk = seq.substring(i, i + KMER_SIZE);

                CompactSerializableInterval csi = new CompactSerializableInterval(rseq.getContigIndex(), i);

                if (!newindex.get(sk).contains(csi)) {
                    throw new IndianaException("Interval '" + csi + "' is not in the index");
                }
            }
        }
    }
}
