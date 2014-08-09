package uk.ac.ox.well.indiana.attic.analyses.LongContigs;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.PrintStream;
import java.util.*;

public class RefRecoveredCDF extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public TreeMap<String, FastaSequenceFile> CONTIGS;

    @Argument(fullName="reference", shortName="r", doc="Reference")
    public TreeMap<String, FastaSequenceFile> REFERENCES;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CortexKmer, Set<String>> refKmers = new HashMap<CortexKmer, Set<String>>();
        Map<String, Integer> totalKmers = new HashMap<String, Integer>();

        log.info("Storing reference kmers...");
        for (String refid : REFERENCES.keySet()) {
            log.info("\t{}", refid);

            FastaSequenceFile ref = REFERENCES.get(refid);

            ReferenceSequence rseq;
            while ((rseq = ref.nextSequence()) != null) {
                String seq = new String(rseq.getBases());

                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                    if (!refKmers.containsKey(kmer)) {
                        refKmers.put(kmer, new TreeSet<String>());
                    }

                    refKmers.get(kmer).add(refid);

                    if (!totalKmers.containsKey(refid)) {
                        totalKmers.put(refid, 0);
                    }
                    totalKmers.put(refid, totalKmers.get(refid) + 1);
                }
            }
        }

        TableWriter tw = new TableWriter(out);

        log.info("Computing % genome recovered per contig...");
        for (String cid : CONTIGS.keySet()) {
            FastaSequenceFile contigs = CONTIGS.get(cid);

            List<ReferenceSequence> rseqs = new ArrayList<ReferenceSequence>();

            log.info("\tloading {} contigs...", cid);
            ReferenceSequence rseq;
            while ((rseq = contigs.nextSequence()) != null) {
                if (rseq.length() > 0) {
                    rseqs.add(rseq);
                }
            }

            log.info("\tsorting {} contigs...", rseqs.size());
            Collections.sort(rseqs, new Comparator<ReferenceSequence>() {
                @Override
                public int compare(ReferenceSequence o1, ReferenceSequence o2) {
                    if (o1.length() == o2.length()) { return 0; }
                    return (o1.length() > o2.length()) ? -1 : 1;
                }
            });

            for (int i = 0; i < 10; i++) {
                log.info("\t{}", rseqs.get(i).length());
            }

            Map<String, Integer> counts = new HashMap<String, Integer>();
            Set<CortexKmer> seenKmers = new HashSet<CortexKmer>();
            int numContigs = 1;

            for (int j = 0; j < rseqs.size(); j++) {
                rseq = rseqs.get(j);
                String seq = new String(rseq.getBases());

                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                    if (refKmers.containsKey(kmer) && !seenKmers.contains(kmer)) {
                        for (String refid : refKmers.get(kmer)) {
                            if (!counts.containsKey(refid)) {
                                counts.put(refid, 1);
                            } else {
                                counts.put(refid, counts.get(refid) + 1);
                            }
                        }

                        seenKmers.add(kmer);
                    }
                }

                Map<String, String> te = new LinkedHashMap<String, String>();
                te.put("num_contigs", String.valueOf(numContigs));
                te.put("algorithm", cid);
                for (String refid : REFERENCES.keySet()) {
                    int seen = counts.containsKey(refid) ? counts.get(refid) : 0;
                    int total = totalKmers.get(refid);

                    double pct = 100.0 * (double) seen / (double) total;
                    te.put(refid, String.format("%.2f", pct));
                }

                tw.addEntry(te);

                numContigs++;
            }
        }
    }
}
