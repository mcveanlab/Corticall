package uk.ac.ox.well.cortexjdk.commands.index.reference;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Set;

public class IndexReference extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference sequence")
    public File REF_FILE;

    @Argument(fullName="source", shortName="s", doc="Link source")
    public String SOURCE;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Argument(fullName="nThreads", shortName="t", doc="Number of threads")
    public Integer NUM_THREADS = 2;

    @Override
    public void execute() {
        log.info("Indexing reference");

        File dbFile = KmerLookup.createIndex(REF_FILE, KMER_SIZE, SOURCE, NUM_THREADS, log);

        try {
            KmerLookup kl = new KmerLookup(REF_FILE);
            ReferenceSequence rseq;
            IndexedFastaSequenceFile fa = new IndexedFastaSequenceFile(REF_FILE);

            while ((rseq = fa.nextSequence()) != null) {
                for (int i = 0; i <= rseq.length() - KMER_SIZE - 10; i++) {
                    String fwd = rseq.getBaseString().substring(i, i + KMER_SIZE);

                    Set<Interval> fwdIntervals = kl.findKmer(fwd);
                    log.info("fwd {} {}", fwd, i);
                    for (Interval it : fwdIntervals) {
                        log.info("  {}", it);
                    }

                    String seq = rseq.getBaseString().substring(i, i + KMER_SIZE + 10);
                    Set<Interval> seqIntervals = kl.findKmer(seq);

                    log.info("seq {} {}", seq, i);
                    for (Interval it : seqIntervals) {
                        log.info("  {}", it);
                    }

                    log.info("");
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        log.info("  index written to {}", dbFile);
    }
}
