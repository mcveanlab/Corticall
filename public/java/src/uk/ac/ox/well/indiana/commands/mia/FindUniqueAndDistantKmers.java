package uk.ac.ox.well.indiana.commands.mia;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

public class FindUniqueAndDistantKmers extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference sequence")
    public FastaSequenceFile REFERENCE;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 21;

    @Output
    public PrintStream out;

    private Map<CortexKmer, String> loadDiagnosticKmers() {
        Map<CortexKmer, String> kmers = new HashMap<CortexKmer, String>();
        Set<CortexKmer> kmersToRemove = new HashSet<CortexKmer>();

        log.info("Loading diagnostic kmers...");
        ReferenceSequence rseq;
        while ((rseq = REFERENCE.nextSequence()) != null) {
            String[] name = rseq.getName().split("\\s+");
            String seq = new String(rseq.getBases());

            log.info("  {}", name[0]);

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                CortexKmer ck = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                if (!kmers.containsKey(ck)) {
                    kmers.put(ck, name[0]);
                } else {
                    kmersToRemove.add(ck);
                }
            }
        }

        log.info("  before: {} kmers", kmers.size());

        for (CortexKmer ck : kmers.keySet()) {
            Collection<byte[]> sks = SequenceUtils.generateSequencesWithEditDistance1(ck.getKmerAsBytes());
            //sks.addAll(SequenceUtils.generateSequencesWithEditDistance1(SequenceUtils.reverseComplement(ck.getKmerAsBytes())));

            boolean removeCk = false;
            for (byte[] sk : sks) {
                CortexKmer ck1 = new CortexKmer(sk);

                if (kmers.containsKey(ck1)) {
                    removeCk = true;
                    kmersToRemove.add(ck1);
                }
            }

            if (removeCk) {
                kmersToRemove.add(ck);
            }
        }

        for (CortexKmer ck : kmersToRemove) {
            kmers.remove(ck);
        }

        log.info("   after: {} kmers", kmers.size());

        return kmers;
    }

    @Override
    public void execute() {
        Map<CortexKmer, String> kmers = loadDiagnosticKmers();

        TableWriter tw = new TableWriter(out);

        for (CortexKmer ck : kmers.keySet()) {
            Map<String, String> te = new LinkedHashMap<String, String>();
            te.put("kmer", ck.getKmerAsString());
            te.put("chr", kmers.get(ck));

            tw.addEntry(te);
        }
    }
}
