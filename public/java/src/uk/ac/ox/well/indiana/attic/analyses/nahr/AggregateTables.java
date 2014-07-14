package uk.ac.ox.well.indiana.attic.analyses.nahr;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.File;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Map;

public class AggregateTables extends Module {
    @Argument(fullName="table", shortName="t", doc="Tables")
    public HashSet<File> TABLES;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        log.info("Aggregating tables...");
        TableWriter tw = new TableWriter(out);

        for (File file : TABLES) {
            log.info("  {}", file.getAbsolutePath());

            if (file.length() > 0) {
                TableReader tr = new TableReader(file);

                for (Map<String, String> te : tr) {
                    tw.addEntry(te);
                }
            }
        }
    }
}
