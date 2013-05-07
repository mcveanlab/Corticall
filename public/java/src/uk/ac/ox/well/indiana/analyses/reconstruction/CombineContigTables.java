package uk.ac.ox.well.indiana.analyses.reconstruction;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.utils.TableReader2;
import uk.ac.ox.well.indiana.utils.io.utils.TableWriter2;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

public class CombineContigTables extends Tool {
    @Argument(fullName="contigTable", shortName="ct", doc="Contig tables")
    public ArrayList<File> CONTIG_TABLES;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableWriter2 tw = new TableWriter2(out);

        for (File contigTable : CONTIG_TABLES) {
            log.info("Processing table '{}'", contigTable.getAbsolutePath());

            TableReader2 tr = new TableReader2(contigTable);

            for (Map<String, String> te : tr) {
                tw.addEntry(te);
            }
        }
    }
}
