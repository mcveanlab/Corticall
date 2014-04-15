package uk.ac.ox.well.indiana.attic.analyses.nahr;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;

import java.io.PrintStream;
import java.util.*;

public class FindChimericContigs extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (in FASTA format)")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="referenceSample", shortName="rs", doc="Reference (in Cortex format)")
    public HashMap<String, CortexGraph> REFERENCES;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, Integer> refIndex = new LinkedHashMap<String, Integer>();

        int index = 0;
        for (String refid : REFERENCES.keySet()) {
            refIndex.put(refid, index);
            index++;
        }

        log.info("Loading reference sequences...");
        Map<CortexKmer, boolean[]> kmerMap = new HashMap<CortexKmer, boolean[]>();
        for (String refid : REFERENCES.keySet()) {
            log.info("\t{}", refid);

            CortexGraph cg = REFERENCES.get(refid);

            for (CortexRecord cr : cg) {
                CortexKmer kmer = cr.getKmer();

                if (!kmerMap.containsKey(kmer)) {
                    kmerMap.put(kmer, new boolean[refIndex.size()]);
                }

                boolean[] mask = kmerMap.get(kmer);
                mask[refIndex.get(refid)] = true;
                kmerMap.put(kmer, mask);
            }
        }

        log.info("\t{} kmers total", kmerMap.size());

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

        log.info("Finding chimeric contigs...");
        int chimeras = 0;
        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            Set<Integer> refIndicesSeen = new HashSet<Integer>();
            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                if (kmerMap.containsKey(kmer)) {
                    boolean[] mask = kmerMap.get(kmer);

                    for (int j = 0; j < mask.length; j++) {
                        if (mask[j]) {
                            refIndicesSeen.add(j);
                            break;
                        }
                    }
                }
            }

            if (refIndicesSeen.size() > 1) {
                out.println(">" + rseq.getName());
                out.println(seq);

                chimeras++;
            }
        }

        log.info("\t{} chimeric contigs", chimeras);
    }
}
