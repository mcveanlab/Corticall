package uk.ac.ox.well.indiana.attic.analyses.kmerSharing;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.performance.PerformanceUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class ExtractRefAndAltKmers extends Module {
    @Argument(fullName="relatedSequences", shortName="rs", doc="Related sequences")
    public ArrayList<File> RELATED_SEQUENCES;

    @Argument(fullName="refKmers", shortName="rk", doc="Enable emission of ref kmers")
    public Boolean EMIT_REF_KMERS = true;

    @Argument(fullName="altKmers", shortName="ak", doc="Enable emission of alt kmers")
    public Boolean EMIT_ALT_KMERS = true;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableWriter tw = new TableWriter(out);

        for (File relatedSequences : RELATED_SEQUENCES) {
            String sample = relatedSequences.getName().replaceAll("relatedSequences.", "").replaceAll(".table", "");

            log.info("Processing sample '{}' ({})", sample, PerformanceUtils.getCompactMemoryUsageStats());

            TableReader tr = new TableReader(relatedSequences);

            for (Map<String, String> e : tr) {
                String kmer = e.get("kmer");
                String supernode = e.get("superNode");
                String gene = e.get("genes");

                for (int i = 0; i <= supernode.length() - kmer.length(); i++) {
                    String fw = SequenceUtils.alphanumericallyLowestOrientation(supernode.substring(i, i + kmer.length()));

                    boolean isRefKmer = fw.equals(kmer);

                    if ((isRefKmer && EMIT_REF_KMERS) || (!isRefKmer && EMIT_ALT_KMERS)) {
                        Map<String, String> entry = new HashMap<String, String>();
                        entry.put("kmer", fw);
                        entry.put("supernodeId", String.valueOf(supernode.hashCode()));
                        entry.put("gene", gene);
                        entry.put(sample, "1");
                    }
                }
            }
        }

        out.println(tw);
    }
}