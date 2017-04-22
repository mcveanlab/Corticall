package uk.ac.ox.well.indiana.commands.playground.caller;

import com.google.common.base.Joiner;
import htsjdk.samtools.*;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.commands.playground.index.KmerIndex;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class AnnotateROIs extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="contam", shortName="c", doc="Contam")
    public CortexGraph CONTAM;

    @Argument(fullName="bam", shortName="b", doc="BAM")
    public File BAM_FILE;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CortexKmer, Map<String, Object>> anns = new HashMap<>();

        KmerIndex ki = new KmerIndex(BAM_FILE, GRAPH.getKmerSize(), false);
        SamReader sr = SamReaderFactory.makeDefault().open(BAM_FILE);

        for (CortexRecord cr : ROI) {
            //Map<String, Object> ann = new HashMap<>();
            //ann.put("")

            List<long[]> chunks = ki.find(cr.getKmerAsBytes());

            List<Byte> quals = new ArrayList<>();

            for (long[] c : chunks) {
                SAMFileSpan sfs = new BAMFileSpan(new Chunk(c[0], c[1]));

                SAMRecordIterator recs = sr.indexing().iterator(sfs);

                while (recs.hasNext()) {
                    SAMRecord rec = recs.next();

                    log.info("    {}: {} {}", cr, new Chunk(c[0], c[1]), rec.getSAMString());

                    for (int i = 0; i <= rec.getReadLength() - GRAPH.getKmerSize(); i++) {
                        CortexKmer ck = new CortexKmer(rec.getReadString().substring(i, i + GRAPH.getKmerSize()));

                        if (cr.getCortexKmer().equals(ck)) {
                            byte minQual = 40;

                            for (int j = i; j < i + GRAPH.getKmerSize(); j++) {
                                byte qual = rec.getBaseQualities()[j];

                                if (qual < minQual) {
                                    minQual = qual;
                                }
                            }

                            quals.add(minQual);
                        }
                    }
                }

                recs.close();
            }

            log.info("{} {}", cr, Joiner.on(", ").join(quals));
        }
    }
}
