package uk.ac.ox.well.indiana.analyses.kmerSharing;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.utils.TableReader;
import uk.ac.ox.well.indiana.utils.io.utils.TableWriter;
import uk.ac.ox.well.indiana.utils.performance.PerformanceUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Map;

public class ExtractRefAndAltKmers extends Tool {
    @Argument(fullName="relatedSequences", shortName="rs", doc="Related sequences")
    public ArrayList<File> RELATED_SEQUENCES;

    @Output
    public PrintStream out;

    @Override
    public int execute() {
        TableWriter tw = new TableWriter();

        tw.addPrimaryKey("kmer", "unknown");
        tw.addColumn("supernodeId", 0);
        tw.addColumn("gene", "unknown");

        for (File relatedSequences : RELATED_SEQUENCES) {
            String sample = relatedSequences.getName().replaceAll("relatedSequences.", "").replaceAll(".table", "");

            log.info("Processing sample '{}' ({})", sample, PerformanceUtils.getCompactMemoryUsageStats());

            tw.addColumn(sample, 0);

            TableReader tr = new TableReader(relatedSequences);

            for (Map<String, String> e : tr) {
                String kmer = e.get("kmer");
                String supernode = e.get("superNode");
                String gene = e.get("genes");

                for (int i = 0; i <= supernode.length() - kmer.length(); i++) {
                    String fw = SequenceUtils.getAlphanumericallyLowestOrientation(supernode.substring(i, i + kmer.length()));

                    tw.set(fw, sample, 1);
                    tw.set(fw, "supernodeId", supernode.hashCode());
                    tw.set(fw, "gene", gene);
                }
            }
        }

        out.println(tw);

        return 0;
    }
}