package uk.ac.ox.well.indiana.analyses.reconstruction;

import uk.ac.ox.well.indiana.tools.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Map;

public class CombineContigTables extends Module {
    @Argument(fullName="contigTable", shortName="ct", doc="Contig tables")
    public ArrayList<File> CONTIG_TABLES;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableWriter tw = new TableWriter(out);

        for (File contigTable : CONTIG_TABLES) {
            log.info("Processing table '{}'", contigTable.getAbsolutePath());

            TableReader tr = new TableReader(contigTable);

            for (Map<String, String> te : tr) {
                tw.addEntry(te);
            }
        }
    }
}
