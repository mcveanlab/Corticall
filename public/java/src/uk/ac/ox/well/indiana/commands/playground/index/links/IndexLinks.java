package uk.ac.ox.well.indiana.commands.playground.index.links;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.mapdb.BTreeMap;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.Serializer;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexBinaryKmer;

import java.io.File;

public class IndexLinks extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference sequence")
    public File REF_FILE;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Override
    public void execute() {
        File dbFile = new File(REF_FILE.getAbsoluteFile() + ".linkdb");

        DB db = DBMaker
                .fileDB(dbFile)
                .fileMmapEnable()
                .make();

        BTreeMap<long[], int[]> kmerIndex = db.treeMap("index")
                .keySerializer(Serializer.LONG_ARRAY)
                .valueSerializer(Serializer.INT_ARRAY)
                .create();

        FastaSequenceFile ref = new FastaSequenceFile(REF_FILE, true);
        ReferenceSequence rseq;
        while ((rseq = ref.nextSequence()) != null) {
            log.info("  {}", rseq.getName());

            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                String sk = seq.substring(i, i + KMER_SIZE);

                if (!sk.contains("N") && !sk.contains(".")) {
                    CortexBinaryKmer cbk = new CortexBinaryKmer(sk.getBytes());

                    int[] l = {rseq.getContigIndex(), i};

                    if (!kmerIndex.containsKey(cbk.getBinaryKmer())) {
                        kmerIndex.put(cbk.getBinaryKmer(), l);
                    } else {
                        int[] l1 = kmerIndex.get(cbk.getBinaryKmer());
                        int[] l2 = new int[l1.length + 2];

                        System.arraycopy(l1, 0, l2, 0, l1.length);
                        System.arraycopy(l, 0, l2, l1.length, l.length);

                        kmerIndex.put(cbk.getBinaryKmer(), l2);
                    }
                }
            }

            db.commit();
        }

        db.close();
    }
}
