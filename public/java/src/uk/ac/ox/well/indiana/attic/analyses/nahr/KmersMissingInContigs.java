package uk.ac.ox.well.indiana.attic.analyses.nahr;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

public class KmersMissingInContigs extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="graph", shortName="g", doc="Graph")
    public HashMap<String, CortexGraph> GRAPHS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, Map<String, Integer>> mkCount = new HashMap<String, Map<String, Integer>>();

        int kmerSize = GRAPHS.values().iterator().next().getKmerSize();

        log.info("Processing contigs...");
        int numContigs = 0;

        TableWriter tw = new TableWriter(out);
        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            if (numContigs % 1000 == 0) {
                log.info("  {} contigs", numContigs);
            }
            numContigs++;

            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                String sk = seq.substring(i, i + kmerSize);
                CortexKmer ck = new CortexKmer(sk);

                for (String label : GRAPHS.keySet()) {
                    CortexGraph cg = GRAPHS.get(label);

                    if (cg.findRecord(ck) == null) {
                        if (!mkCount.containsKey(rseq.getName())) {
                            mkCount.put(rseq.getName(), new HashMap<String, Integer>());
                        }

                        if (!mkCount.get(rseq.getName()).containsKey(label)) {
                            mkCount.get(rseq.getName()).put(label, 0);
                        }

                        mkCount.get(rseq.getName()).put(label, 1 + mkCount.get(rseq.getName()).get(label));
                    }
                }
            }

            Map<String, String> te = new HashMap<String, String>();
            for (String label : GRAPHS.keySet()) {
                te.put("contigName", rseq.getName());
                te.put(label, "0");

                if (mkCount.containsKey(rseq.getName()) && mkCount.get(rseq.getName()).containsKey(label)) {
                    te.put(label, String.valueOf(mkCount.get(rseq.getName()).get(label)));
                }
            }

            tw.addEntry(te);
        }
    }
}
