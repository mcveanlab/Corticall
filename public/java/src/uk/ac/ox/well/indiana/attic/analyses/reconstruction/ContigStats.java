package uk.ac.ox.well.indiana.attic.analyses.reconstruction;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class ContigStats extends Module {
    @Argument(fullName="contigTable", shortName="ct", doc="Contig tables")
    public ArrayList<File> CONTIG_TABLES;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableWriter tw = new TableWriter(out);

        for (File table : CONTIG_TABLES) {
            Map<String, Collection<String>> contigsPerSample = new HashMap<String, Collection<String>>();

            TableReader tr = new TableReader(table);

            for (Map<String, String> te : tr) {
                String sample = te.get("sample");
                String contig = te.get("contig");

                if (!contigsPerSample.containsKey(sample)) {
                    contigsPerSample.put(sample, new ArrayList<String>());
                }

                contigsPerSample.get(sample).add(contig);
            }

            for (String sample : contigsPerSample.keySet()) {
                Map<String, String> te = new LinkedHashMap<String, String>();
                int minLength = SequenceUtils.minLength(contigsPerSample.get(sample));
                int maxLength = SequenceUtils.maxLength(contigsPerSample.get(sample));
                int n50Length = SequenceUtils.computeN50Length(contigsPerSample.get(sample));
                int n50Value = SequenceUtils.computeN50Value(contigsPerSample.get(sample));

                te.put("sample", sample);
                te.put("numContigs", String.valueOf(contigsPerSample.get(sample).size()));
                te.put("minLength", String.valueOf(maxLength));
                te.put("maxLength", String.valueOf(maxLength));
                te.put("n50Length", String.valueOf(n50Length));
                te.put("n50Value", String.valueOf(n50Value));

                tw.addEntry(te);
            }
        }
    }
}
