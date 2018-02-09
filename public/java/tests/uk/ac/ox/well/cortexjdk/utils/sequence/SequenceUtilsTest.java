package uk.ac.ox.well.cortexjdk.utils.sequence;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.gff.GFF3;
import uk.ac.ox.well.cortexjdk.utils.io.gff.GFF3Record;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class SequenceUtilsTest {
    @Test
    public void complementTest() {
        byte[] trials = new byte[] { 'A', 'C', 'G', 'T', 'N', '.', 'a', 'c', 'g', 't' };
        byte[] exp    = new byte[] { 'T', 'G', 'C', 'A', 'N', '.', 't', 'g', 'c', 'a' };

        for (int i = 0; i < trials.length; i++) {
            Assert.assertEquals(SequenceUtils.complement(trials[i]), exp[i]);
        }
    }

    @Test
    public void reverseComplementTest() {
        byte[] sequence   = "TACTGACTTTTCTCGCTATTCGTATGCATG".getBytes();
        byte[] expectedRC = "CATGCATACGAATAGCGAGAAAAGTCAGTA".getBytes();

        byte[] computedRC = SequenceUtils.reverseComplement(sequence);

        Assert.assertEquals(expectedRC, computedRC);
    }

    @Test
    public void reverseComplementTestWithNs() {
        byte[] sequence   = "NACTGACTTTTCTCGCTATTCGTATGCATG".getBytes();
        byte[] expectedRC = "CATGCATACGAATAGCGAGAAAAGTCAGTN".getBytes();

        byte[] computedRC = SequenceUtils.reverseComplement(sequence);

        Assert.assertEquals(expectedRC, computedRC);
    }

    @Test
    public void reverseComplementWithLowercaseBases() {
        byte[] sequence   = "NACTGACTTTTCTCGCTATTCGTATGCATg".getBytes();
        byte[] expectedRC = "cATGCATACGAATAGCGAGAAAAGTCAGTN".getBytes();

        byte[] computedRC = SequenceUtils.reverseComplement(sequence);

        Assert.assertEquals(expectedRC, computedRC);
    }

    @Test
    public void alphanumericallyLowestOrientation() {
        for (int i = 0; i < 10000; i++) {
            for (Integer k : Arrays.asList(21, 31, 41, 51)) {
                String fw = new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(k));
                String rc = SequenceUtils.reverseComplement(fw);

                byte[] expected = fw.compareTo(rc) < 0 ? fw.getBytes() : rc.getBytes();
                byte[] computed = SequenceUtils.alphanumericallyLowestOrientation(fw.getBytes());

                Assert.assertEquals(expected, computed);
            }
        }
    }

    @Test
    public void computeN50() {
        List<String> sequences = new ArrayList<>();
        sequences.add("AAGCTTA");
        sequences.add("TTGA");
        sequences.add("AAC");
        sequences.add("TT");
        sequences.add("AA");
        sequences.add("C");
        sequences.add("G");

        int n50 = SequenceUtils.computeN50Value(sequences);

        Assert.assertEquals(4, n50);

        List<CanonicalKmer> csequences = new ArrayList<>();
        csequences.add(new CanonicalKmer("AAGCTTA"));
        csequences.add(new CanonicalKmer("TTGA"));
        csequences.add(new CanonicalKmer("AAC"));
        csequences.add(new CanonicalKmer("TT"));
        csequences.add(new CanonicalKmer("AA"));
        csequences.add(new CanonicalKmer("C"));
        csequences.add(new CanonicalKmer("G"));

        int cn50 = SequenceUtils.computeN50Value(csequences);

        Assert.assertEquals(4, cn50);
    }

    @Test(enabled = false)
    public void testTranslateCodingSequence() {
        String cds1 = "ATGGTGACTGGTAGTGGTGGTGAGGATAAGTATAAAAGTGCCAAAAATGCCAAGGAACTT";
        String tr1  = "MVTGSGGEDKYKSAKNAKEL";

        Assert.assertEquals(tr1, SequenceUtils.translateCodingSequence(cds1));
    }
}
