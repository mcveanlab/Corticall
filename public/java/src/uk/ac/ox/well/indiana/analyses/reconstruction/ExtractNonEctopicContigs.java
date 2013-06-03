package uk.ac.ox.well.indiana.analyses.reconstruction;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class ExtractNonEctopicContigs extends Tool {
    @Argument(fullName="contigTable", shortName="ct", doc="Contig table")
    public File CONTIG_TABLE;

    @Argument(fullName="sampleName", shortName="sn", doc="Sample to extract")
    public HashSet<String> SAMPLES;

    @Argument(fullName="geneName", shortName="gn", doc="Gene names")
    public HashSet<String> GENES;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableReader tr = new TableReader(CONTIG_TABLE);

        int id = 0;
        for (Map<String, String> te : tr) {
            String sample = te.get("sample");

            Set<String> genes = new HashSet<String>(Arrays.asList(te.get("genes").split(",")));

            if (SAMPLES.contains(sample) && genes.size() == 1) {
                for (String geneOfInterest : GENES) {
                    if (genes.contains(geneOfInterest)) {
                        String contig = te.get("contig");
                        String contigName = ">contig_" + contig.hashCode();

                        out.println(contigName);
                        out.println(te.get("contig"));

                        break;
                    }
                }
            }

            id++;
        }
    }
}
