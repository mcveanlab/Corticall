package uk.ac.ox.well.indiana.attic.analyses.nahr.old;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;

import java.io.PrintStream;
import java.util.*;

public class KmerOverlapWithReference extends Module {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public ArrayList<CortexGraph> CORTEX_GRAPHS;

    @Argument(fullName="reference", shortName="R", doc="Reference sequences")
    public HashMap<String, IndexedFastaSequenceFile> REFERENCES;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int kmerSize = CORTEX_GRAPHS.get(0).getKmerSize();

        Map<String, Integer> refIndex = new LinkedHashMap<String, Integer>();

        int index = 0;
        for (String refid : REFERENCES.keySet()) {
            refIndex.put(refid, index);
            index++;
        }

        log.info("Loading reference sequences...");
        Map<CortexKmer, boolean[]> kmerMap = new HashMap<CortexKmer, boolean[]>();
        for (String refid : REFERENCES.keySet()) {
            log.info("  {}", refid);

            IndexedFastaSequenceFile ref = REFERENCES.get(refid);

            ReferenceSequence rseq;
            while ((rseq = ref.nextSequence()) != null) {
                String seq = new String(rseq.getBases());

                for (int i = 0; i <= seq.length() - kmerSize; i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + kmerSize));

                    if (!kmerMap.containsKey(kmer)) {
                        kmerMap.put(kmer, new boolean[refIndex.size()]);
                    }

                    boolean[] mask = kmerMap.get(kmer);
                    mask[refIndex.get(refid)] = true;
                }
            }
        }

        log.info("  {} kmers total", kmerMap.size());

        log.info("Removing kmers found in multiple references...");
        Set<CortexKmer> kmersToRemove = new HashSet<CortexKmer>();
        for (CortexKmer kmer : kmerMap.keySet()) {
            boolean[] mask = kmerMap.get(kmer);

            int presenceCount = 0;
            for (boolean value : mask) {
                if (value) {
                    presenceCount++;
                }
            }

            if (presenceCount > 1) {
                kmersToRemove.add(kmer);
            }
        }

        for (CortexKmer kmer : kmersToRemove) {
            kmerMap.remove(kmer);
        }

        log.info("\t{} kmers to remove, {} kmers remaining", kmersToRemove.size(), kmerMap.size());

        kmersToRemove.clear();

        Map<String, Integer> uniqueKmers = new HashMap<String, Integer>();
        for (CortexKmer kmer : kmerMap.keySet()) {
            boolean[] mask = kmerMap.get(kmer);

            for (String refid : refIndex.keySet()) {
                int refi = refIndex.get(refid);

                if (mask[refi]) {
                    if (!uniqueKmers.containsKey(refid)) {
                        uniqueKmers.put(refid, 1);
                    } else {
                        uniqueKmers.put(refid, uniqueKmers.get(refid) + 1);
                    }
                }
            }
        }

        log.info("{}", uniqueKmers);

        log.info("Processing graphs...");

        out.println("sample" + "\t" + Joiner.on("\t").join(refIndex.keySet()));

        for (CortexGraph cg : CORTEX_GRAPHS) {
            Map<String, Integer> kmerCounts = new HashMap<String, Integer>();

            for (String refName : refIndex.keySet()) {
                kmerCounts.put(refName, 0);
            }

            for (CortexRecord cr : cg) {
                CortexKmer kmer = cr.getKmer();

                if (kmerMap.containsKey(kmer)) {
                    boolean[] mask = kmerMap.get(kmer);

                    for (String refid : refIndex.keySet()) {
                        int refi = refIndex.get(refid);

                        if (mask[refi]) {
                            kmerCounts.put(refid, kmerCounts.get(refid) + 1);
                        }
                    }
                }
            }

            List<String> counts = new ArrayList<String>();
            for (String refName : refIndex.keySet()) {
                counts.add(String.valueOf((float) kmerCounts.get(refName) / (float) uniqueKmers.get(refName)));
            }

            out.println(cg.getColor(0).getSampleName() + "\t" + Joiner.on("\t").join(counts));
        }
    }
}
