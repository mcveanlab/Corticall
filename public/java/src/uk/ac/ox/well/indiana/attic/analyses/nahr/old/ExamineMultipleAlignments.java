package uk.ac.ox.well.indiana.attic.analyses.nahr.old;

import com.google.common.base.Joiner;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;

import java.io.PrintStream;
import java.util.*;

public class ExamineMultipleAlignments extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public ArrayList<SAMFileReader> BAMS;

    @Argument(fullName="referenceSample", shortName="rs", doc="Reference (in Cortex format)")
    public HashMap<String, CortexGraph> REFERENCES;

    //@Argument(fullName="reference", shortName="r", doc="Reference (in FASTA format)")
    //public HashMap<String, IndexedFastaSequenceFile> FASTAS;

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

        log.info("Loading reads...");

        Map<String, Set<SAMRecord>> reads = new HashMap<String, Set<SAMRecord>>();

        for (SAMFileReader b : BAMS) {
            log.info("  {}", b);

            for (SAMRecord r : b) {
                String readName = r.getReadName();
                if (!reads.containsKey(readName)) {
                    reads.put(readName, new HashSet<SAMRecord>());
                }

                reads.get(readName).add(r);
            }
        }

        int kmerSize = REFERENCES.values().iterator().next().getKmerSize();

        for (String readName : reads.keySet()) {
            String seq = reads.get(readName).iterator().next().getReadString();
            StringBuilder painting = new StringBuilder();
            painting.setLength(seq.length());

            for (int i = 0; i < seq.length(); i++) {
                painting.setCharAt(i, '-');
            }

            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + kmerSize));

                if (kmerMap.containsKey(kmer)) {
                    boolean[] mask = kmerMap.get(kmer);

                    for (int j = 0; j < mask.length; j++) {
                        if (mask[j]) {
                            for (int k = 0; k < kmerSize; k++) {
                                if (painting.charAt(i + k) == '-') {
                                    painting.setCharAt(i + k, String.valueOf(j).charAt(0));
                                }
                            }

                            break;
                        }
                    }
                }
            }

            Set<CigarElement> ces = new HashSet<CigarElement>();

            out.println(readName + ":");

            out.println("  seq: " + seq);
            out.println("  pai: " + new String(painting));
            out.println("  pa2: " + (new String(painting)).replaceAll("-", ""));

            for (SAMRecord r : reads.get(readName)) {
                out.print("  " + r.getSAMString());

                for (CigarElement ce : r.getCigar().getCigarElements()) {
                    ces.add(ce);
                }
            }

            List<String> cesStrings = new ArrayList<String>();
            for (CigarElement ce : ces) {
                cesStrings.add(ce.getOperator() + " " + ce.getLength());
            }

            out.println("  " + Joiner.on("; ").join(cesStrings));

            out.println();
        }
    }
}
