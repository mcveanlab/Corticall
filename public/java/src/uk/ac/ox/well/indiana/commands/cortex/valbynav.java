package uk.ac.ox.well.indiana.commands.cortex;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Description;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

@Description(text="Validates contigs by navigating the available Cortex graphs")
public class valbynav extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (FASTA)")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="contigName", shortName="n", doc="Select specific contigs to process by name", required=false)
    public HashSet<String> CONTIG_NAMES;

    @Argument(fullName="cortexGraphs", shortName="g", doc="Sorted cortex graphs")
    public HashSet<CortexGraph> CORTEX_GRAPHS;

    @Argument(fullName="reverseComplement", shortName="rc", doc="Reverse-complement the selected contigs")
    public Boolean REVERSE_COMPLEMENT = false;

    @Output
    public PrintStream out;

    private Set<ReferenceSequence> loadSequences() {
        Set<ReferenceSequence> sequences = new LinkedHashSet<ReferenceSequence>();

        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            String[] name = rseq.getName().split("\\s+");

            if (CONTIG_NAMES == null || CONTIG_NAMES.isEmpty() || CONTIG_NAMES.contains(name[0])) {
                sequences.add(rseq);
            }
        }

        return sequences;
    }

    @Override
    public void execute() {
        log.info("Loading contigs...");
        Set<ReferenceSequence> sequences = loadSequences();
        log.info("  loaded {} contigs", sequences.size());

        for (ReferenceSequence rseq : sequences) {
            String seq = new String(rseq.getBases());

            if (REVERSE_COMPLEMENT) {
                seq = SequenceUtils.reverseComplement(seq);
            }

            Map<Integer, String> output = new TreeMap<Integer, String>();
            String skLast = null;

            for (CortexGraph cg : CORTEX_GRAPHS) {
                int kmerSize = cg.getKmerSize();

                StringBuilder statusCodes = new StringBuilder();

                for (int i = 0; i <= seq.length() - kmerSize; i++) {
                    String sk = seq.substring(i, i + kmerSize);
                    CortexKmer ck = new CortexKmer(sk);

                    CortexRecord cr = cg.findRecord(ck);

                    String code = " ";

                    if (cr == null) { code = "."; }
                    if (cr != null && skLast != null) {
                        String kmerNext = CortexUtils.getNextKmer(cg, skLast);
                        String kmerPrev = CortexUtils.getPrevKmer(cg, sk);

                        if (kmerNext != null && kmerPrev != null && (!kmerNext.equals(sk) || !kmerPrev.equals(skLast))) {
                            code = "x";
                        }
                    }

                    statusCodes.append(code);

                    skLast = sk;
                }

                output.put(kmerSize, statusCodes.toString());
            }

            log.info("");
            log.info("  name: {}", rseq.getName());
            log.info("contig: {}", seq);
            for (int kmerSize : output.keySet()) {
                log.info("   k{}: {}", kmerSize, output.get(kmerSize));
            }
        }
    }
}
