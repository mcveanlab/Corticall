package uk.ac.ox.well.indiana.commands.playground.links;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.PrintStream;

public class FindKmersFromContigs extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="dirty", shortName="d", doc="Dirty graph")
    public CortexGraph GRAPH;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int kmersUsedTotal = 0;
        int contigsTotal = 0;

        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            String contigName = rseq.getName().split("\\s+")[0];

            String seq = rseq.getBaseString();

            log.info("contig: {}", contigName);

            int contigs = 0;
            int counter = 0;

            for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                String sk = seq.substring(i, i + GRAPH.getKmerSize());
                CortexKmer ck = new CortexKmer(sk);
                CortexRecord cr = GRAPH.findRecord(ck);

                kmersUsedTotal++;

                if (cr == null) {
                    if (counter > 0) {
                        log.info("{}", sk);

                        for (int j = i - 5; j < i; j++) {
                            String sk0 = seq.substring(j, j + GRAPH.getKmerSize());
                            CortexKmer ck0 = new CortexKmer(sk0);
                            CortexRecord cr0 = GRAPH.findRecord(ck0);

                            log.info("{} {} {} {}", contigName, j, sk0, cr0);
                        }

                        int limit = seq.length() - GRAPH.getKmerSize();
                        boolean limitReset = false;
                        for (int j = i; j < limit; j++) {
                            String sk0 = seq.substring(j, j + GRAPH.getKmerSize());
                            CortexKmer ck0 = new CortexKmer(sk0);
                            CortexRecord cr0 = GRAPH.findRecord(ck0);

                            log.info("{} {} {} {}", contigName, j, sk0, cr0);

                            if (cr0 != null && !limitReset) {
                                limitReset = true;
                                limit = j + 5;
                            }

                            if (limitReset && cr0 == null) {
                                limit = j + 5;
                            }
                        }

                        if (limitReset) {
                            i += limit;
                        }

                        //log.info("counter: {}", counter);
                        counter = 0;
                        contigs++;

                        log.info("contigs: {}", contigs);
                        log.info("");
                    }
                } else {
                    counter++;
                }
            }

            log.info("  counter={} contigs={}", counter, contigs);
            contigsTotal += contigs;
        }

        log.info("kmers={} contigs={}", kmersUsedTotal, contigsTotal);
    }
}
