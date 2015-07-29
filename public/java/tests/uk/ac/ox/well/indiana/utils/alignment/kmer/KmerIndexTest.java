package uk.ac.ox.well.indiana.utils.alignment.kmer;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class KmerIndexTest {
    @Test
    public void testLookup() {
        File refFile = new File("testdata/twochrs.fasta");
        FastaSequenceFile ref = new FastaSequenceFile(refFile, true);
        ReferenceSequence rseq;
        int kmerSize = 47;

        Map<String, Set<Interval>> locations = new HashMap<String, Set<Interval>>();

        while ((rseq = ref.nextSequence()) != null) {
            String name = rseq.getName().split("\\s+")[0];
            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                String sk = seq.substring(i, i + kmerSize);
                Interval interval = new Interval(name, i, i);

                if (!locations.containsKey(sk)) {
                    locations.put(sk, new HashSet<Interval>());
                }

                locations.get(sk).add(interval);
            }
        }

        KmerIndex ki = new KmerIndex(refFile, kmerSize);

        for (String sk : locations.keySet()) {
            KmerIndex.KmerIndexRecord kir = ki.lookup(sk);

            Set<Interval> kilocs = kir.getLocations();
            Set<Interval> locs = locations.get(sk);

            Assert.assertEquals(kilocs, locs);
        }
    }
}
