package uk.ac.ox.well.indiana.commands.attic;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.io.PrintStream;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

/**
 * Created by kiran on 04/08/2016.
 */
public class EmitKmers extends Module {
    @Argument(fullName="seq", shortName="s", doc="Kmers")
    public File SEQ;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 30;

    @Argument(fullName="limit", shortName="l", doc="Limit")
    public Integer LIMIT = 10;

    @Output
    public PrintStream out;

    private class ValueComparator implements Comparator<String> {
        Map<String, Integer> base;

        public ValueComparator(Map<String, Integer> base) {
            this.base = base;
        }

        // Note: this comparator imposes orderings that are inconsistent with
        // equals.
        public int compare(String a, String b) {
            if (base.get(a) >= base.get(b)) {
                return -1;
            } else {
                return 1;
            } // returning 0 would merge keys
        }
    }

    @Override
    public void execute() {
        //ReferenceSequence rseq;
        StringBuilder sb = new StringBuilder();
        LineReader lr = new LineReader(SEQ);
        while (lr.hasNext()) {
            sb.append(lr.getNextRecord());
        }

        String rseq = sb.toString();

        Map<String, Integer> unsorted = new HashMap<>();
        ValueComparator bvc = new ValueComparator(unsorted);
        Map<String, Integer> sorted = new TreeMap<>(bvc);

        for (int i = 0; i <= rseq.length() - KMER_SIZE; i++) {
            String seq = rseq.substring(i, i + KMER_SIZE);

            ContainerUtils.increment(unsorted, seq);
        }

        sorted.putAll(unsorted);

        log.info("kmers: {}", unsorted.size());

        int i = 0;
        for (String sk : sorted.keySet()) {
            //log.info("{} {} {}", i, sk, unsorted.get(sk));

            out.printf("%s,%d\n", sk, unsorted.get(sk));

            i++;

            if (i >= LIMIT) { break; }
        }
    }
}
