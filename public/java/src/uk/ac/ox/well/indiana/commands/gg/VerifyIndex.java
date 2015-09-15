package uk.ac.ox.well.indiana.commands.gg;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.mapdb.BTreeKeySerializer;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.Fun;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.NavigableSet;

public class VerifyIndex extends Module {
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
                .readOnly()
                .make();

        NavigableSet<Object[]> kmerIndex = db.treeSet("index" + KMER_SIZE);

        try {
            log.info("Verifying kmer locations...");

            IndexedFastaSequenceFile fa = new IndexedFastaSequenceFile(REF_FILE);

            int index = 0;
            for (Object[] l : kmerIndex) {
                if (index % (kmerIndex.size() / 10) == 0) {
                    log.info("  {}/{} kmers", index, kmerIndex.size());
                }

                String sk = (String) l[0];
                String chr = (String) l[1];
                int pos = (Integer) l[2];

                String kmer = new String(fa.getSubsequenceAt(chr, pos + 1, pos + KMER_SIZE).getBases());

                if (!sk.equals(kmer)) {
                    throw new IndianaException("Indexed kmer '" + sk + "' has position " + chr + ":" + pos + ", but kmer found at that position doesn't match ('" + kmer + "')");
                }

                index++;
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}
