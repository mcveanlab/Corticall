package uk.ac.ox.well.indiana.commands.gg;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTables;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class NovelStretchesInGraph extends Module {
    @Argument(fullName = "graphClean", shortName = "c", doc = "Graph (clean)")
    public CortexGraph CLEAN;

    @Argument(fullName = "graphDirty", shortName = "d", doc = "Graph (dirty)", required = false)
    public CortexGraph DIRTY;

    @Argument(fullName = "novelGraph", shortName = "n", doc = "Graph of novel kmers")
    public CortexGraph NOVEL;

    @Argument(fullName="skipToKmer", shortName="s", doc="Skip processing to given kmer", required=false)
    public String KMER;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CortexKmer, Boolean> novelKmers = new HashMap<CortexKmer, Boolean>();
        for (CortexRecord cr : NOVEL) {
            novelKmers.put(new CortexKmer(cr.getKmerAsString()), true);
        }

        int totalNovelKmersUsed = 0;
        int stretchNum = 1;

        Set<CortexKmer> novelKmersToVisit = new HashSet<CortexKmer>();

        if (KMER != null) {
            novelKmersToVisit.add(new CortexKmer(KMER));
        } else {
            novelKmersToVisit.addAll(novelKmers.keySet());
        }

        TableWriter tw = new TableWriter(out);

        log.info("Discovering novel kmer stretches in graph...");
        for (CortexKmer novelKmer : novelKmersToVisit) {
            if (novelKmers.get(novelKmer)) {
                // Walk the graph left and right of novelKmer and extract a novel stretch
                String stretch;
                if (DIRTY == null) {
                    stretch = CortexUtils.getSeededStretch(CLEAN, novelKmer.getKmerAsString(), 0, true);
                } else {
                    stretch = CortexUtils.getSeededStretch(CLEAN, DIRTY, novelKmer.getKmerAsString(), 0, true);
                }

                int novelKmersUsed = 0;

                for (int i = 0; i <= stretch.length() - CLEAN.getKmerSize(); i++) {
                    String sk = stretch.substring(i, i + CLEAN.getKmerSize());
                    CortexKmer ck = new CortexKmer(sk);

                    if (novelKmers.containsKey(ck)) {
                        novelKmersUsed++;
                        totalNovelKmersUsed++;
                        novelKmers.put(ck, false);
                    }
                }

                log.info("  stretch {}:", stretchNum);
                log.info("    novel kmer:             {}", novelKmer);
                log.info("    length (bp):            {}", stretch.length());
                log.info("    novel kmers:            {}", novelKmersUsed);
                log.info("    total novel kmers used: {}", totalNovelKmersUsed);
                log.info("    sequence:");
                log.info("    - 0: {}", SequenceUtils.truncate(stretch, 100));

                Map<String, String> te = new LinkedHashMap<String, String>();
                te.put("stretchNum", String.valueOf(stretchNum));
                te.put("novelKmer", novelKmer.getKmerAsString());
                te.put("length", String.valueOf(stretch.length()));
                te.put("novelKmers", String.valueOf(novelKmersUsed));
                te.put("cumNovelKmers", String.valueOf(totalNovelKmersUsed));
                te.put("totalNovelKmers", String.valueOf(novelKmersToVisit.size()));
                te.put("stretch", stretch);

                tw.addEntry(te);

                stretchNum++;
            }
        }
    }
}
