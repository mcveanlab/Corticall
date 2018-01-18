package uk.ac.ox.well.cortexjdk.utils.alignment.kmer;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
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
    private IndexedReference kl;

    private Map<String, Set<Interval>> expectedIntervals = new HashMap<>();
    private Map<Interval, String> expectedKmersFwd = new HashMap<>();
    private Map<Interval, String> expectedKmersRev = new HashMap<>();
    private int[] expectedKmerSizes = { 31, 47 };
    private String expectedSource = "test";

    public KmerLookupTest() {
        File testFa = new File("tests/two_short_contigs.fa");
        File testFai = new File(testFa.getAbsolutePath() + ".fai");
        File testDict = new File(testFa.getAbsolutePath() + ".dict");

        File testAmb = new File(testFa.getAbsolutePath() + ".amb");
        File testAnn = new File(testFa.getAbsolutePath() + ".ann");
        File testBwt = new File(testFa.getAbsolutePath() + ".bwt");
        File testPac = new File(testFa.getAbsolutePath() + ".pac");
        File testSa  = new File(testFa.getAbsolutePath() + ".sa");

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

            File tempAmb = new File(tempFa.getAbsolutePath() + ".amb");
            File tempAnn = new File(tempFa.getAbsolutePath() + ".ann");
            File tempBwt = new File(tempFa.getAbsolutePath() + ".bwt");
            File tempPac = new File(tempFa.getAbsolutePath() + ".pac");
            File tempSa  = new File(tempFa.getAbsolutePath() + ".sa");

            FileUtils.copyFile(testFa, tempFa);
            FileUtils.copyFile(testFai, tempFai);
            FileUtils.copyFile(testDict, tempDict);

            FileUtils.copyFile(testAmb, tempAmb);
            FileUtils.copyFile(testAnn, tempAnn);
            FileUtils.copyFile(testBwt, tempBwt);
            FileUtils.copyFile(testPac, tempPac);
            FileUtils.copyFile(testSa, tempSa);

            tempFa.deleteOnExit();
            tempFai.deleteOnExit();
            tempDict.deleteOnExit();

            IndexedReference.createIndex(tempFa, expectedSource).deleteOnExit();

            kl = new IndexedReference(tempFa);
        } catch (IOException e) {
            throw new CortexJDKException("Could not write KmerLookupTest files to temp directory");
        }
    }

    @Test
    public void testSource() {
        Assert.assertTrue(kl.getSources().contains("test"));
    }

    @Test
    public void findKmerBySequence() {
        for (String sk : expectedIntervals.keySet()) {
            System.out.println(sk + " " + expectedIntervals.get(sk) + " " + kl.find(sk));
            Assert.assertEquals(expectedIntervals.get(sk), kl.find(sk));
        }
    }

    @Test
    public void findKmerByIntervalFwd() {
        for (Interval it : expectedKmersFwd.keySet()) {
            Assert.assertEquals(expectedKmersFwd.get(it), kl.find(it));
        }
    }

    @Test
    public void findKmerByIntervalRev() {
        for (Interval it : expectedKmersRev.keySet()) {
            Assert.assertEquals(expectedKmersRev.get(it), kl.find(it));
        }
    }

    @Test
    public void findMissingKmerReturnsEmptyList() {
        Assert.assertEquals(0, kl.find("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT").size());
    }
}
