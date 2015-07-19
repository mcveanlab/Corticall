package uk.ac.ox.well.indiana.commands.gg;

import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class PrintNovelStretches extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="novelGraph", shortName="n", doc="Graph of novel kmers")
    public CortexGraph NOVEL;

    @Argument(fullName="novelKmerMap", shortName="m", doc="Novel kmer map")
    public File NOVEL_KMER_MAP;

    @Output
    public PrintStream out;

    private class VariantInfo {
        public String variantId;
        public String vclass;
        public String vchr;
        public int vstart;
        public int vstop;

        @Override
        public String toString() {
            return "VariantInfo{" +
                    "variantId='" + variantId + '\'' +
                    ", vclass='" + vclass + '\'' +
                    ", vchr='" + vchr + '\'' +
                    ", vstart=" + vstart +
                    ", vstop=" + vstop +
                    '}';
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            VariantInfo that = (VariantInfo) o;

            if (vstart != that.vstart) return false;
            if (vstop != that.vstop) return false;
            if (variantId != null ? !variantId.equals(that.variantId) : that.variantId != null) return false;
            if (vclass != null ? !vclass.equals(that.vclass) : that.vclass != null) return false;
            return !(vchr != null ? !vchr.equals(that.vchr) : that.vchr != null);

        }

        @Override
        public int hashCode() {
            int result = variantId != null ? variantId.hashCode() : 0;
            result = 31 * result + (vclass != null ? vclass.hashCode() : 0);
            result = 31 * result + (vchr != null ? vchr.hashCode() : 0);
            result = 31 * result + vstart;
            result = 31 * result + vstop;
            return result;
        }
    }

    private Map<CortexKmer, VariantInfo> loadNovelKmerMap() {
        Map<CortexKmer, VariantInfo> vis = new HashMap<CortexKmer, VariantInfo>();

        TableReader tr = new TableReader(NOVEL_KMER_MAP);

        for (Map<String, String> te : tr) {
            VariantInfo vi = new VariantInfo();

            if (!te.get("vstart").equals("NA")) {
                vi.variantId = te.get("variantId");
                vi.vclass = te.get("vclass");
                vi.vchr = te.get("vchr");
                vi.vstart = Integer.valueOf(te.get("vstart"));
                vi.vstop = Integer.valueOf(te.get("vstop"));

                CortexKmer kmer = new CortexKmer(te.get("kmer"));

                vis.put(kmer, vi);
            }
        }

        return vis;
    }

    @Override
    public void execute() {
        Map<CortexKmer, VariantInfo> vis = loadNovelKmerMap();

        Map<CortexKmer, Boolean> novelKmers = new HashMap<CortexKmer, Boolean>();
        for (CortexRecord cr : NOVEL) {
            novelKmers.put(new CortexKmer(cr.getKmerAsString()), true);
        }

        int totalNovelKmersUsed = 0;
        int stretchNum = 0;
        for (CortexKmer novelKmer : novelKmers.keySet()) {
            if (novelKmers.get(novelKmer)) {
                int novelKmersUsed = 0;

                String stretch = CortexUtils.getSeededStretch(GRAPH, novelKmer.getKmerAsString(), 0);
                Set<VariantInfo> vs = new HashSet<VariantInfo>();

                for (int i = 0; i <= stretch.length() - novelKmer.length(); i++) {
                    CortexKmer ck = new CortexKmer(stretch.substring(i, i + novelKmer.length()));

                    if (novelKmers.containsKey(ck) && novelKmers.get(ck)) {
                        totalNovelKmersUsed++;
                        novelKmersUsed++;
                        novelKmers.put(ck, false);
                    }

                    if (vis.containsKey(ck)) {
                        vs.add(vis.get(ck));
                    }
                }

                for (VariantInfo vi : vs) {
                    Set<CortexKmer> entriesToRemove = new HashSet<CortexKmer>();
                    for (CortexKmer ck : vis.keySet()) {
                        if (vi.equals(vis.get(ck))) {
                            entriesToRemove.add(ck);
                        }
                    }

                    for (CortexKmer ck : entriesToRemove) {
                        vis.remove(ck);
                    }
                }

                log.info("Found {} bp stretch, used {}/{} novel kmers, {}/{} total",
                    stretch.length(),
                    novelKmersUsed, novelKmers.size(),
                    totalNovelKmersUsed, novelKmers.size()
                );

                for (VariantInfo vi : vs) {
                    log.info("\t{}", vi);
                }

                out.println(">stretch_" + stretchNum + ".length_" + stretch.length());
                out.println(stretch);

                stretchNum++;
            }
        }

        Set<VariantInfo> remainingVariants = new HashSet<VariantInfo>(vis.values());
        log.info("Remaining variants: {}", remainingVariants.size());
        for (VariantInfo vi : remainingVariants) {
            log.info("  {}", vi);
        }
    }
}
