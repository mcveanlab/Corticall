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

public class AnnotateContigs extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (in FASTA format)")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="referenceSample", shortName="rs", doc="Reference (in Cortex format)")
    public HashMap<String, CortexGraph> REFERENCES;

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

        /*
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
        */

        int kmerSize = kmerMap.keySet().iterator().next().length();

        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            StringBuilder annotation = new StringBuilder();

            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + kmerSize));

                if (!kmerMap.containsKey(kmer)) {
                    annotation.append(".");
                } else {
                    boolean[] mask = kmerMap.get(kmer);
                    int presenceIndex = 0;
                    int presenceCount = 0;

                    for (int j = 0; j < mask.length; j++) {
                        if (mask[j]) {
                            presenceCount++;
                            presenceIndex = j;
                        }
                    }

                    if (presenceCount == 1) {
                        annotation.append(String.valueOf(presenceIndex));
                    } else {
                        annotation.append("B");
                    }
                }
            }

            out.println(rseq.getName() + "\t" + seq + "\t" + annotation.toString());
        }
    }
}
