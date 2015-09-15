package uk.ac.ox.well.indiana.commands.gg;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.mapdb.BTreeKeySerializer;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Map;
import java.util.NavigableSet;
import java.util.Set;
import java.util.concurrent.ConcurrentNavigableMap;

public class IndexReference extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference sequence")
    public File REF_FILE;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        File dbFile = new File(REF_FILE.getAbsoluteFile() + ".kmerdb");

        DB db = DBMaker.fileDB(dbFile)
                .closeOnJvmShutdown()
                .transactionDisable()
                .fileMmapEnable()
                .cacheSize(1000000)
                .make();

        NavigableSet<Object[]> kmerIndex = db.treeSetCreate("index" + KMER_SIZE)
                .serializer(BTreeKeySerializer.ARRAY3)
                .make();

        FastaSequenceFile ref = new FastaSequenceFile(REF_FILE, true);
        ReferenceSequence rseq;
        while ((rseq = ref.nextSequence()) != null) {
            log.info("  {}", rseq.getName());

            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                String sk = seq.substring(i, i + KMER_SIZE);

                //kmerIndex.add(new Object[]{sk, rseq.getContigIndex(), i});
                kmerIndex.add(new Object[]{sk, rseq.getName().split("\\s+")[0], i});
            }

            db.commit();
        }

        db.close();
    }
}
