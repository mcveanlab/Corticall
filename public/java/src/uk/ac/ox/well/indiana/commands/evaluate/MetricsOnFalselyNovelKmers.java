package uk.ac.ox.well.indiana.commands.evaluate;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.PrintStream;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Set;

public class MetricsOnFalselyNovelKmers extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Output
    public PrintStream out;

    private class StretchInfo {
        int stretchNum;
        boolean isOrphaned = true;
        String stretch;
        int novelKmers = 0;

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            StretchInfo that = (StretchInfo) o;

            if (stretchNum != that.stretchNum) return false;
            if (isOrphaned != that.isOrphaned) return false;
            if (novelKmers != that.novelKmers) return false;
            return !(stretch != null ? !stretch.equals(that.stretch) : that.stretch != null);

        }

        @Override
        public int hashCode() {
            int result = stretchNum;
            result = 31 * result + (isOrphaned ? 1 : 0);
            result = 31 * result + (stretch != null ? stretch.hashCode() : 0);
            result = 31 * result + novelKmers;
            return result;
        }
    }

    @Override
    public void execute() {
        Set<CortexKmer> usedKmers = new HashSet<CortexKmer>();
        Set<CortexKmer> orphanKmers = new HashSet<CortexKmer>();

        Set<StretchInfo> sis = new LinkedHashSet<StretchInfo>();
        int stretchNum = 0, numOrphaned = 0, numOrphanKmers = 0;

        for (CortexRecord cr : GRAPH) {
            if (cr.getCoverage(0) > 0) {
                if (cr.getCoverage(1) > 0) {
                    //out.println(num + " " + cr);
                } else if (!usedKmers.contains(cr.getCortexKmer())) {
                    StretchInfo si = new StretchInfo();
                    si.stretchNum = stretchNum;
                    si.stretch = CortexUtils.getSeededStretch(GRAPH, cr.getKmerAsString(), 2, true);

                    log.info("seq {} {}:", stretchNum, cr.getKmerAsString());
                    log.info("  {} {}", si.stretch.length(), si.stretch);

                    for (int i = 0; i <= si.stretch.length() - GRAPH.getKmerSize(); i++) {
                        String sk = si.stretch.substring(i, i + GRAPH.getKmerSize());
                        CortexKmer ck = new CortexKmer(sk);
                        CortexRecord rec = GRAPH.findRecord(ck);

                        log.info("  {} {} {}", i, sk, rec);

                        usedKmers.add(ck);

                        if (rec.getCoverage(0) > 0) {
                            si.novelKmers++;
                        }

                        if (rec.getCoverage(3) > 0 || rec.getCoverage(4) > 0) {
                            si.isOrphaned = false;
                        }
                    }

                    log.info("  isOrphaned: {}, novelKmers: {}, length: {}", si.isOrphaned, si.novelKmers, si.stretch.length());

                    if (si.isOrphaned) {
                        numOrphaned++;
                        numOrphanKmers += si.novelKmers;

                        for (int i = 0; i <= si.stretch.length() - GRAPH.getKmerSize(); i++) {
                            String sk = si.stretch.substring(i, i + GRAPH.getKmerSize());
                            CortexKmer ck = new CortexKmer(sk);
                            CortexRecord rec = GRAPH.findRecord(ck);

                            if (rec.getCoverage(0) > 0) {
                                orphanKmers.add(ck);
                            }
                        }
                    }

                    sis.add(si);

                    stretchNum++;
                }
            }
        }

        log.info("     numOrphaned: {}", numOrphaned);
        log.info("numOrphanedKmers: {}", orphanKmers.size());
    }
}
