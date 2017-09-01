package uk.ac.ox.well.cortexjdk.commands.index.reference;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;

import java.io.File;

public class IndexReference extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference sequence")
    public File REF_FILE;

    @Argument(fullName="source", shortName="s", doc="Link source")
    public String SOURCE;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Override
    public void execute() {
        log.info("Indexing reference");

        File dbFile = KmerLookup.createIndex(REF_FILE, SOURCE, KMER_SIZE, log);

        log.info("  index written to {}", dbFile);
    }
}
