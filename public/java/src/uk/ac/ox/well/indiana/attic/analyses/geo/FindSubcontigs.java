package uk.ac.ox.well.indiana.attic.analyses.geo;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

public class FindSubcontigs extends Module {
    @Argument(fullName="contigs1", shortName="c1", doc="Contigs from first sample")
    public FastaSequenceFile CONTIGS1;

    @Argument(fullName="contigs2", shortName="c2", doc="Contigs from second sample")
    public FastaSequenceFile CONTIGS2;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 21;

    @Output
    public PrintStream out;

    private Map<ReferenceSequence, Set<CortexKmer>> loadContigs(FastaSequenceFile contigs) {
        Map<ReferenceSequence, Set<CortexKmer>> contigsAndKmers = new HashMap<ReferenceSequence, Set<CortexKmer>>();

        ReferenceSequence rseq;
        while ((rseq = contigs.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            contigsAndKmers.put(rseq, new HashSet<CortexKmer>());

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                contigsAndKmers.get(rseq).add(kmer);
            }
        }

        return contigsAndKmers;
    }

    private boolean hasKmerOverlap(Set<CortexKmer> kmers1, Set<CortexKmer> kmers2) {
        return !Collections.disjoint(kmers1, kmers2);
    }

    @Override
    public void execute() {
        Map<ReferenceSequence, Set<CortexKmer>> contigsAndKmers1 = loadContigs(CONTIGS1);
        Map<ReferenceSequence, Set<CortexKmer>> contigsAndKmers2 = loadContigs(CONTIGS2);

        //Map<String, Map<String, String>> subContigs = new HashMap<String, Map<String, String>>();

        for (ReferenceSequence rseq1 : contigsAndKmers1.keySet()) {
            //log.info("{}: {}", rseq1.getName(), new String(rseq1.getBases()));

            Map<String, String> subContigNames = new HashMap<String, String>();
            Map<String, String> subContigs = new TreeMap<String, String>(new Comparator<String>() {
                @Override
                public int compare(String o1, String o2) {
                    if (o1.length() == o2.length()) { return 0; }
                    return o1.length() < o2.length() ? 1 : -1;
                }
            });
            String contig1 = new String(rseq1.getBases());

            for (ReferenceSequence rseq2 : contigsAndKmers2.keySet()) {
                if (hasKmerOverlap(contigsAndKmers1.get(rseq1), contigsAndKmers2.get(rseq2))) {
                    String fw2 = new String(rseq2.getBases());
                    String rc2 = SequenceUtils.reverseComplement(fw2);

                    String lcsFw = SequenceUtils.longestCommonSubstring(contig1, fw2);
                    String lcsRc = SequenceUtils.longestCommonSubstring(contig1, rc2);
                    String lcs = lcsFw.length() > lcsRc.length() ? lcsFw : lcsRc;
                    String contig2 = lcsFw.length() > lcsRc.length() ? fw2 : rc2;

                    subContigs.put(lcs, contig2);
                    subContigNames.put(lcs, rseq2.getName());
                }
            }

            if (!subContigs.isEmpty()) {
                Map.Entry<String, String> entry = subContigs.entrySet().iterator().next();
                String lcs = entry.getKey();
                String contig2 = entry.getValue();

                out.println(rseq1.getName() + " " +
                            contig1 + " " +
                            subContigNames.get(lcs) + " " +
                            contig2 + " " +
                            lcs + " " +
                            lcs.length()
                );
            }
        }
    }
}
