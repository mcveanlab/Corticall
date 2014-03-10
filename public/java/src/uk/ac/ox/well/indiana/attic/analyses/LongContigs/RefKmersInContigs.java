package uk.ac.ox.well.indiana.attic.analyses.LongContigs;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

public class RefKmersInContigs extends Module {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex binary")
    public HashMap<String, CortexGraph> CORTEX_GRAPHS;

    @Argument(fullName="contigs", shortName="c", doc="Contigs FASTA")
    public HashMap<String, FastaSequenceFile> CONTIGS;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        log.info("Storing reference kmers...");

        Map<CortexKmer, Map<String, Boolean>> refKmersSeen = new HashMap<CortexKmer, Map<String,Boolean>>();
        Map<String, Integer> refKmersTotal = new HashMap<String, Integer>();
        for (String refid : CORTEX_GRAPHS.keySet()) {
            log.info("  {}...", refid);

            CortexGraph ref = CORTEX_GRAPHS.get(refid);
            for (CortexRecord cr : ref) {
                CortexKmer kmer = cr.getKmer();

                if (!refKmersSeen.containsKey(kmer)) {
                    refKmersSeen.put(kmer, new HashMap<String, Boolean>());
                }

                if (!refKmersSeen.get(kmer).containsKey(refid)) {
                    refKmersSeen.get(kmer).put(refid, false);

                    if (!refKmersTotal.containsKey(refid)) {
                        refKmersTotal.put(refid, 1);
                    } else {
                        refKmersTotal.put(refid, refKmersTotal.get(refid) + 1);
                    }
                }
            }
        }

        Map<String, Integer> refKmersUnique = new HashMap<String, Integer>();
        for (CortexKmer kmer : refKmersSeen.keySet()) {
            if (refKmersSeen.get(kmer).size() == 1) {
                for (String refid : refKmersSeen.get(kmer).keySet()) {
                    if (!refKmersUnique.containsKey(refid)) {
                        refKmersUnique.put(refid, 1);
                    } else {
                        refKmersUnique.put(refid, refKmersUnique.get(refid) + 1);
                    }
                }
            }
        }

        for (String refid : refKmersUnique.keySet()) {
            log.info("  {} : {} total kmers, {} unique kmers", refid, refKmersTotal.get(refid), refKmersUnique.get(refid));
        }

        log.info("Examining contigs...");

        TableWriter tw = new TableWriter(out);

        for (String cid : CONTIGS.keySet()) {
            log.info("  {}...", cid);

            FastaSequenceFile contig = CONTIGS.get(cid);

            Set<String> contigs = new HashSet<String>();

            ReferenceSequence rseq;
            while ((rseq = contig.nextSequence()) != null) {
                String seq = new String(rseq.getBases());

                contigs.add(seq);

                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                    if (refKmersSeen.containsKey(kmer)) {
                        for (String refid : refKmersSeen.get(kmer).keySet()) {
                            refKmersSeen.get(kmer).put(refid, true);
                        }
                    }
                }
            }

            Map<String, Integer> seen = new TreeMap<String, Integer>();
            Map<String, Integer> unique = new TreeMap<String, Integer>();
            for (CortexKmer kmer : refKmersSeen.keySet()) {
                for (String refid : refKmersSeen.get(kmer).keySet()) {
                    if (refKmersSeen.get(kmer).get(refid)) {
                        if (!seen.containsKey(refid)) {
                            seen.put(refid, 1);
                        } else {
                            seen.put(refid, seen.get(refid) + 1);
                        }

                        if (refKmersSeen.get(kmer).size() == 1) {
                            if (!unique.containsKey(refid)) {
                                unique.put(refid, 1);
                            } else {
                                unique.put(refid, unique.get(refid) + 1);
                            }
                        }

                        refKmersSeen.get(kmer).put(refid, false);
                    }
                }
            }

            for (String refid : seen.keySet()) {
                Map<String, String> te = new LinkedHashMap<String, String>();

                te.put("id", cid);
                te.put("num_contigs", String.valueOf(contigs.size()));
                te.put("n50", String.valueOf(SequenceUtils.computeN50Value(contigs)));
                te.put("ref", refid);
                te.put("total", String.valueOf(refKmersTotal.get(refid)));
                te.put("unique", String.valueOf(refKmersUnique.get(refid)));
                te.put("seen_total", String.valueOf(seen.get(refid)));

                if (unique.containsKey(refid)) {
                    te.put("seen_unique", String.valueOf(unique.get(refid)));
                } else {
                    te.put("seen_unique", String.valueOf(0.0));
                }

                te.put("pct_seen_total", String.valueOf(100.0 * seen.get(refid) / (refKmersTotal.get(refid))));

                if (unique.containsKey(refid)) {
                    te.put("pct_seen_unique", String.valueOf(100.0 * unique.get(refid) / (refKmersUnique.get(refid))));
                } else {
                    te.put("pct_seen_unique", String.valueOf(0.0));
                }

                tw.addEntry(te);
            }
        }
    }
}
