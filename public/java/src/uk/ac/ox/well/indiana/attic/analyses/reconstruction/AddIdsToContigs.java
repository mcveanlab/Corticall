package uk.ac.ox.well.indiana.attic.analyses.reconstruction;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.File;
import java.io.PrintStream;
import java.util.Map;

public class AddIdsToContigs extends Module {
    @Argument(fullName="contigTable", shortName="ct", doc="Contig table")
    public File CONTIG_TABLE;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableReader tr = new TableReader(CONTIG_TABLE);
        TableWriter tw = new TableWriter(out);

        for (Map<String, String> te : tr) {
            String contig = te.get("contig");
            te.put("contigId", String.valueOf(contig.hashCode()));

            tw.addEntry(te);
        }
    }
}
