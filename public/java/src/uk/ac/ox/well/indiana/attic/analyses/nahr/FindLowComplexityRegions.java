package uk.ac.ox.well.indiana.attic.analyses.nahr;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

public class FindLowComplexityRegions extends Module {
    @Argument(fullName="fasta", shortName="f", doc="Fasta")
    public IndexedFastaSequenceFile FASTA;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 41;

    @Argument(fullName="padding", shortName="p", doc="Padding")
    public Integer PADDING = 33;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, Set<Interval>> kmers = new HashMap<String, Set<Interval>>();

        Set<String> allowableContigs = new HashSet<String>();
        allowableContigs.add("Pf3D7_02_v3");
        allowableContigs.add("Pf3D7_01_v3");

        Map<String, Integer> chrLength = new HashMap<String, Integer>();

        log.info("Processing genome...");
        ReferenceSequence rseq;
        while ((rseq = FASTA.nextSequence()) != null) {
            if (allowableContigs.contains(rseq.getName())) {
                log.info("  {}", rseq.getName());

                String seq = new String(rseq.getBases());

                chrLength.put(rseq.getName(), seq.length());

                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    String kmer = seq.substring(i, i + KMER_SIZE);
                    if (!kmers.containsKey(kmer)) {
                        kmers.put(kmer, new HashSet<Interval>());
                    }

                    if (kmers.get(kmer).size() <= 3) {
                        String[] name = rseq.getName().split("\\s+");
                        Interval interval = new Interval(name[0], i + 1, i + 1 + KMER_SIZE);

                        kmers.get(kmer).add(interval);
                    }
                }
            }
        }

        for (String kmer : kmers.keySet()) {
            if (kmers.get(kmer).size() == 2) {
                Iterator<Interval> regions = kmers.get(kmer).iterator();
                Interval ia = regions.next();
                Interval ib = regions.next();

                if (!ia.getSequence().equals(ib.getSequence())) {
                    Interval ia1 = new Interval(ia.getSequence(), ia.getStart() - PADDING, ia.getEnd() + PADDING);
                    Interval ib1 = new Interval(ib.getSequence(), ib.getStart() - PADDING, ib.getEnd() + PADDING);

                    if (ia1.getStart() > 0 && ib1.getStart() > 0 && ia1.getEnd() < chrLength.get(ia1.getSequence()) && ib1.getEnd() < chrLength.get(ib1.getSequence())) {
                        ReferenceSequence ra = FASTA.getSubsequenceAt(ia1.getSequence(), ia1.getStart(), ia1.getEnd());
                        ReferenceSequence rb = FASTA.getSubsequenceAt(ib1.getSequence(), ib1.getStart(), ib1.getEnd());

                        String sa = new String(ra.getBases());
                        String sb = new String(rb.getBases());

                        String sa0 = sa.substring(0, PADDING);
                        String sb0 = sb.substring(0, PADDING);

                        String sa1 = sa.substring(PADDING, sa.length() - PADDING);
                        String sb1 = sb.substring(PADDING, sb.length() - PADDING);

                        String sa2 = sa.substring(sa.length() - PADDING, sa.length() - 1);
                        String sb2 = sb.substring(sb.length() - PADDING, sb.length() - 1);

                        if (!sa0.equals(sb0) && sa1.equals(sb1) && !sa2.equals(sb2) && SequenceUtils.editDistance(sa2.getBytes(), sb2.getBytes()) > 10) {
                            out.println(kmer + "\t" + Joiner.on("\t").join(ia1, ib1));
                        }
                    }
                }
            }
        }
    }
}
