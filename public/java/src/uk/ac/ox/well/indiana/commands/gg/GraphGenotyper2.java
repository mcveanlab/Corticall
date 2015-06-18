package uk.ac.ox.well.indiana.commands.gg;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class GraphGenotyper2 extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="novelKmerGraph", shortName="n", doc="Novel kmer graph")
    public CortexGraph NOVEL_KMER_GRAPH;

    @Argument(fullName="aggressivelyAggregateNovelKmers", shortName="a", doc="Be aggressive in novel kmer aggregation")
    public Boolean AGGRESSIVE = false;

    @Output
    public PrintStream out;

    public String walkUntilWeSeeParentalJunctions(CortexGraph cg, String novelKmer, int color, Map<CortexKmer, Boolean> usableNovelKmers) {
        String tk = novelKmer;
        StringBuilder stretchBuilder = new StringBuilder(tk);

        do {
            Set<String> pks = CortexUtils.getPrevKmers(cg, tk, color);

            tk = null;

            if (pks.size() == 1) {
                tk = pks.iterator().next();
            } else if (pks.size() > 1 && AGGRESSIVE) {
                for (String pk0 : pks) {
                    if (usableNovelKmers.containsKey(new CortexKmer(pk0))) {
                        tk = pk0;
                    }
                }
            }

            if (tk != null) {
                stretchBuilder.insert(0, String.valueOf(tk.charAt(0)));
            }
        } while (tk != null);

        tk = novelKmer;

        do {
            Set<String> nks = CortexUtils.getNextKmers(cg, tk, color);

            tk = null;

            if (nks.size() == 1) {
                tk = nks.iterator().next();
            } else if (nks.size() > 1 && AGGRESSIVE) {
                for (String nk0 : nks) {
                    if (usableNovelKmers.containsKey(new CortexKmer(nk0))) {
                        tk = nk0;
                    }
                }
            }

            if (tk != null) {
                stretchBuilder.append(String.valueOf(tk.charAt(tk.length() - 1)));
            }
        } while (tk != null);

        return stretchBuilder.toString();
    }

    @Override
    public void execute() {
        Map<CortexKmer, Boolean> usableNovelKmers = new HashMap<CortexKmer, Boolean>();

        for (CortexRecord cr : NOVEL_KMER_GRAPH) {
            CortexKmer ck = cr.getCortexKmer();

            usableNovelKmers.put(ck, true);
        }

        int numStretches = 0;
        for (CortexKmer ck : usableNovelKmers.keySet()) {
            if (usableNovelKmers.get(ck)) {
                String stretch = walkUntilWeSeeParentalJunctions(GRAPH, ck.getKmerAsString(), 0, usableNovelKmers);

                String fk = stretch.substring(0, NOVEL_KMER_GRAPH.getKmerSize());
                String lk = stretch.substring(stretch.length() - NOVEL_KMER_GRAPH.getKmerSize(), stretch.length());

                log.info("{} {}: {}", numStretches, ck, stretch.length());

                numStretches++;

                for (int i = 0; i <= stretch.length() - NOVEL_KMER_GRAPH.getKmerSize(); i++) {
                    CortexKmer nk = new CortexKmer(stretch.substring(i, i + NOVEL_KMER_GRAPH.getKmerSize()));

                    log.info("  {}/{} ({}): {}", i, stretch.length() - NOVEL_KMER_GRAPH.getKmerSize(), usableNovelKmers.containsKey(nk) ? "n" : "k", GRAPH.findRecord(nk));

                    if (usableNovelKmers.containsKey(nk)) {
                        usableNovelKmers.put(nk, false);
                    }
                }
            }
        }
    }
}
