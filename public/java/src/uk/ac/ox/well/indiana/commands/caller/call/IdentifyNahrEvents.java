package uk.ac.ox.well.indiana.commands.caller.call;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by kiran on 22/06/2017.
 */
public class IdentifyNahrEvents extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="refs", shortName="R", doc="References")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Argument(fullName="sequences", shortName="s", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="validated", shortName="v", doc="Validated NAHR event")
    public FastaSequenceFile NAHR;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CortexKmer, Set<ReferenceSequence>> expectedNovelKmers = new HashMap<>();
        ReferenceSequence rseq;
        while ((rseq = NAHR.nextSequence()) != null) {
            String seq = rseq.getBaseString();
            for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                CortexKmer ck = new CortexKmer(seq.substring(i, i + GRAPH.getKmerSize()));
                CortexRecord rr = ROI.findRecord(ck);

                if (rr != null) {
                    expectedNovelKmers.put(rr.getCortexKmer(), new HashSet<>());
                }
            }
        }

        while ((rseq = CONTIGS.nextSequence()) != null) {
            String seq = rseq.getBaseString();
            for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                CortexKmer ck = new CortexKmer(seq.substring(i, i + GRAPH.getKmerSize()));

                if (expectedNovelKmers.containsKey(ck)) {
                    ContainerUtils.add(expectedNovelKmers, ck, rseq);
                }
            }
        }

        for (CortexKmer ck : expectedNovelKmers.keySet()) {
            log.info("ck={} rseqs={}", ck, expectedNovelKmers.get(ck));
        }

    }
}
