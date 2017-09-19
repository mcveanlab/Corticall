package uk.ac.ox.well.cortexjdk.utils.alignment.kmer;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * Created by kiran on 31/08/2017.
 */
public class KmerLookupTest {
    private File tempFa;
    private KmerLookup kl;

    private Map<String, Set<Interval>> expectedIntervals = new HashMap<>();
    private Map<Interval, String> expectedKmersFwd = new HashMap<>();
    private Map<Interval, String> expectedKmersRev = new HashMap<>();
    private int[] expectedKmerSizes = { 3, 5 };
    private String expectedSource = "test";

    public KmerLookupTest() {
        File testFa = new File("tests/two_short_contigs.fa");
        File testFai = new File(testFa.getAbsolutePath() + ".fai");
        File testDict = new File(testFa.getAbsolutePath() + ".dict");

        try {
            IndexedFastaSequenceFile ref = new IndexedFastaSequenceFile(testFa);

            ReferenceSequence rseq;
            while ((rseq = ref.nextSequence()) != null) {
                String seq = rseq.getBaseString();

                for (int kmerSize : expectedKmerSizes) {
                    for (int i = 0; i <= seq.length() - kmerSize; i++) {
                        String kmerFwd = seq.substring(i, i + kmerSize);
                        String kmerRev = SequenceUtils.reverseComplement(kmerFwd);

                        Interval iFwd = new Interval(rseq.getName(), i + 1, i + kmerSize, false, null);
                        Interval iRev = new Interval(rseq.getName(), i + 1, i + kmerSize, true, null);

                        if (!expectedIntervals.containsKey(kmerFwd)) { expectedIntervals.put(kmerFwd, new HashSet<>()); }
                        expectedIntervals.get(kmerFwd).add(iFwd);
                        expectedKmersFwd.put(iFwd, kmerFwd);

                        if (!expectedIntervals.containsKey(kmerRev)) { expectedIntervals.put(kmerRev, new HashSet<>()); }
                        expectedIntervals.get(kmerRev).add(iRev);
                        expectedKmersRev.put(iRev, kmerRev);
                    }
                }
            }
        } catch (FileNotFoundException e) {
            throw new CortexJDKException("KmerLookupTest data not found.");
        }

        try {
            tempFa = File.createTempFile("KmerLookupTest", ".fa");

            File tempFai = new File(tempFa.getAbsolutePath() + ".fai");
            File tempDict = new File(tempFa.getAbsolutePath() + ".dict");

            FileUtils.copyFile(testFa, tempFa);
            FileUtils.copyFile(testFai, tempFai);
            FileUtils.copyFile(testDict, tempDict);

            tempFa.deleteOnExit();
            tempFai.deleteOnExit();
            tempDict.deleteOnExit();

            for (int kmerSize : expectedKmerSizes) {
                KmerLookup.createIndex(tempFa, kmerSize, expectedSource, 1, null).deleteOnExit();
            }

            kl = new KmerLookup(tempFa);
        } catch (IOException e) {
            throw new CortexJDKException("Could not write KmerLookupTest files to temp directory");
        }
    }

    @Test
    public void testKmerSizes() {
        Assert.assertEquals(kl.getKmerSizes(), new HashSet<>(Arrays.asList(3, 5)));
    }

    @Test
    public void testSource() {
        Assert.assertEquals(kl.getSource(), "test");
    }

    @Test
    public void findKmerBySequence() {
        for (String sk : expectedIntervals.keySet()) {
            Assert.assertEquals(expectedIntervals.get(sk), kl.findKmer(sk));
        }
    }

    @Test
    public void findKmerByIntervalFwd() {
        for (Interval it : expectedKmersFwd.keySet()) {
            Assert.assertEquals(expectedKmersFwd.get(it), kl.findKmer(it));
        }
    }

    @Test
    public void findKmerByIntervalRev() {
        for (Interval it : expectedKmersRev.keySet()) {
            Assert.assertEquals(expectedKmersRev.get(it), kl.findKmer(it));
        }
    }

    @Test
    public void findMissingKmerReturnsEmptyList() {
        Assert.assertEquals(0, kl.findKmer("GTACC").size());
    }

    @Test(expectedExceptions = CortexJDKException.class)
    public void findKmerSizeWithNoCorrespondingKmerdbThrowsException() {
        kl.findKmer("TCTGCATATGA");
    }

    @Test(expectedExceptions = CortexJDKException.class)
    public void findIntervalOnNonExistentContigThrowsException() {
        kl.findKmer(new Interval("3", 1, 2));
    }
}
