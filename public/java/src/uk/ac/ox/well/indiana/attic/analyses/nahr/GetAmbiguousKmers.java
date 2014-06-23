package uk.ac.ox.well.indiana.attic.analyses.nahr;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.Indiana;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

public class GetAmbiguousKmers extends Module {
    @Argument(fullName="fasta", shortName="f", doc="FASTA")
    public IndexedFastaSequenceFile FASTA;

    @Argument(fullName="mask", shortName="m", doc="Mask")
    public File MASK;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Output
    public PrintStream out;

    @Output(fullName="statsOut", shortName="so", doc="Stats out")
    public PrintStream sout;

    @Override
    public void execute() {
        Map<String, Integer> kmerCategories = new HashMap<String, Integer>();

        LineReader lr = new LineReader(MASK);
        String line;
        int lineNum = 0;
        while ((line = lr.getNextRecord()) != null) {
            if (lineNum >= 3) {
                String chr, type;
                int start, stop;

                try {
                    String[] fields = line.split("\\s+");
                    chr = fields[5];
                    start = Integer.valueOf(fields[6]);
                    stop = Integer.valueOf(fields[7]);
                    type = fields[11];
                } catch (NumberFormatException e) {
                    log.info("{}: {}", lineNum, line);

                    throw new IndianaException("Failed to parse line" + e);
                }

                log.info("chr={} start={} end={} type={}", chr, start, stop, type);

                ReferenceSequence rseq = FASTA.getSubsequenceAt(chr, start, stop);
                String seq = new String(rseq.getBases());

                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + KMER_SIZE));

                    String name = String.format("%s:%d-%d.%d %s", chr, start, stop, i, type);
                    out.println(">" + name);
                    out.println(kmer.getKmerAsString());

                    if (!kmerCategories.containsKey(type)) {
                        kmerCategories.put(type, 0);
                    }
                    kmerCategories.put(type, kmerCategories.get(type) + 1);
                }
            }

            lineNum++;
        }

        for (String key : kmerCategories.keySet()) {
            sout.printf("%s\t%d\n", key, kmerCategories.get(key));
        }
    }
}
