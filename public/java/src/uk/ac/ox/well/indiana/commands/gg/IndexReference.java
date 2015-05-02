package uk.ac.ox.well.indiana.commands.gg;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.kmerindex.KmerIndex;

import java.io.*;

public class IndexReference extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public File REFERENCE;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Override
    public void execute() {
        log.info("Indexing reference...");
        int numRecordsWritten = KmerIndex.writeIndex(REFERENCE, KMER_SIZE);

        log.info("  wrote {} records", numRecordsWritten);
    }
}
