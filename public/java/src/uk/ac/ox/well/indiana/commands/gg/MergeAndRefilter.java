package uk.ac.ox.well.indiana.commands.gg;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTable;
import uk.ac.ox.well.indiana.utils.containers.DataTables;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class MergeAndRefilter extends Module {
    @Argument(fullName="clean", shortName="c", doc="Graph (clean)")
    public CortexGraph CLEAN;

    @Argument(fullName="dirty", shortName="d", doc="Graph (dirty)")
    public CortexGraph DIRTY;

    @Argument(fullName="novelKmers", shortName="n", doc="Novel kmer graph")
    public CortexGraph NOVEL_KMERS;

    @Argument(fullName="calls", shortName="calls", doc="Calls")
    public File CALLS;

    @Argument(fullName="alignments", shortName="a", doc="Alignments")
    public SAMFileReader ALIGNMENTS;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Argument(fullName="covThreshold", shortName="t", doc="Coverage threshold")
    public Integer COV_THRESHOLD = 5;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<CortexKmer> novelKmers = new HashSet<CortexKmer>();
        for (CortexRecord cr : NOVEL_KMERS) {
            novelKmers.add(cr.getCortexKmer());
        }

        Map<CortexKmer, SAMRecord> srs = new HashMap<CortexKmer, SAMRecord>();

        for (SAMRecord sr : ALIGNMENTS) {
            String sk = sr.getReadName().split("\\.")[0];
            CortexKmer ck = new CortexKmer(sk);

            srs.put(ck, sr);
        }

        DataTables dts = new DataTables(CALLS);

        DataTable dt = dts.getTable("variantCalls");

        for (String pk : dt.getPrimaryKeys()) {
            Map<String, Object> de = dt.get(pk);

            String sk = (String) de.get("novelKmer");
            CortexKmer ck = new CortexKmer(sk);

            String childStretch = (String) de.get("childStretch");

            if (childStretch != null) {
                boolean pass = true;

                for (int i = 0; i <= childStretch.length() - KMER_SIZE - 1; i++) {
                    String csk = childStretch.substring(i,     i + KMER_SIZE);
                    String nsk = childStretch.substring(i + 1, i + KMER_SIZE + 1);

                    String nextEdge = nsk.substring(nsk.length() - 1, nsk.length());

                    CortexKmer cck = new CortexKmer(csk);

                    CortexRecord ccr = CLEAN.findRecord(cck);
                    CortexRecord cdr = DIRTY.findRecord(cck);

                    boolean cleanParentHasEdge = ccr.getOutEdgesAsStrings(1).contains(nextEdge) || ccr.getOutEdgesAsStrings(2).contains(nextEdge);
                    boolean dirtyParentHasEdge = cdr.getOutEdgesAsStrings(1).contains(nextEdge) || cdr.getOutEdgesAsStrings(2).contains(nextEdge);

                    if (ccr.getCoverage(0) > COV_THRESHOLD && CortexUtils.isNovelKmer(ccr, 0) && !novelKmers.contains(cck) && !cleanParentHasEdge && dirtyParentHasEdge) {
                        pass = false;
                    }
                }

                if (!pass) {
                    de.put("filter", "PARENTS_OVERCLEANED");
                }

                if (srs.get(ck) != null) {
                    SAMRecord sr = srs.get(ck);

                    String locus = String.format("%s:%d", sr.getReferenceName(), sr.getAlignmentStart());

                    de.put("locus", locus);
                }

                dt.set(pk, de);

                log.info("{} {} {} {}", ck, pass, de, srs.get(ck));
            }
        }

        dts.write(out);
    }
}
