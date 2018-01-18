package uk.ac.ox.well.cortexjdk.commands.index.reference;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;

import java.io.File;
import java.util.ArrayList;

public class IndexReference extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference sequence")
    public File REF_FILE;

    @Argument(fullName="source", shortName="s", doc="Link source")
    public ArrayList<String> SOURCES;

    @Override
    public void execute() {
        log.info("Indexing reference");

        File dbFile = IndexedReference.createIndex(REF_FILE, SOURCES.toArray(new String[0]));

        log.info("  wrote {} sources to {}", SOURCES.size(), dbFile.getAbsolutePath());
    }
}
