package uk.ac.ox.well.indiana.attic.analyses.LongContigs;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.PrintStream;
import java.util.*;

public class ClassifyContigs extends Module {
    @Argument(fullName="fasta", shortName="f", doc="FASTA file")
    public FastaSequenceFile FASTA;

    @Argument(fullName="reference", shortName="r", doc="Reference FASTA files")
    public HashMap<String, IndexedFastaSequenceFile> REFERENCES;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        log.info("Loading kmers from references...");
        Map<CortexKmer, Set<String>> rkmers = new HashMap<CortexKmer, Set<String>>();
        Map<String, Integer> seqids = new HashMap<String, Integer>();

        ReferenceSequence rseq;
        int idnum = 1;
        for (String id : REFERENCES.keySet()) {
            seqids.put(id, idnum);
            idnum++;

            IndexedFastaSequenceFile ref = REFERENCES.get(id);

            while ((rseq = ref.nextSequence()) != null) {
                //log.info("    {}", rseq.getName());

                String seq = new String(rseq.getBases());

                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                    if (!rkmers.containsKey(kmer)) {
                        rkmers.put(kmer, new TreeSet<String>());
                    }

                    rkmers.get(kmer).add(id);
                }
            }
        }

        log.info("{}", seqids);

        log.info("Classifying contigs...");
        int processed = 0;

        TableWriter tw = new TableWriter(out);
        while ((rseq = FASTA.nextSequence()) != null) {
            if (processed % 10000 == 0) {
                //log.info("  Processed {} contigs", processed);
            }
            processed++;

            String seq = new String(rseq.getBases());

            if (seq.length() > 0) {
                Map<String, String> entry = new LinkedHashMap<String, String>();
                entry.put("name", rseq.getName());

                Map<String, Integer> classifications = new TreeMap<String, Integer>();

                StringBuilder cb = new StringBuilder();

                boolean saw1 = false;
                boolean saw2 = false;
                boolean sawMix = false;

                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                    if (rkmers.containsKey(kmer)) {
                        if (rkmers.get(kmer).size() == 1) {
                            for (String id : rkmers.get(kmer)) {
                                if (!classifications.containsKey(id)) {
                                    classifications.put(id, 0);
                                }

                                classifications.put(id, classifications.get(id) + 1);

                                cb.append(seqids.get(id));

                                if (seqids.get(id) == 1) { saw1 = true; }
                                if (seqids.get(id) == 2) { saw2 = true; }
                            }
                        } else {
                            cb.append("?");

                            sawMix = true;
                        }
                    } else {
                        cb.append(" ");
                    }
                }

                if (saw1 && saw2) {
                    log.info("{}", rseq.getName());
                    log.info("{}", seq);
                    log.info("{}", cb.toString());
                    log.info("");
                }

                for (String id : REFERENCES.keySet()) {
                    if (classifications.containsKey(id)) {
                        entry.put(id, String.valueOf(classifications.get(id)));
                    } else {
                        entry.put(id, String.valueOf(0));
                    }
                }

                entry.put("contig", seq);

                tw.addEntry(entry);
            }
        }
    }
}
