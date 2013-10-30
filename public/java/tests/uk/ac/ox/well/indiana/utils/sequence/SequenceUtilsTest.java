package uk.ac.ox.well.indiana.utils.sequence;

import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;

import java.util.ArrayList;
import java.util.List;

public class SequenceUtilsTest {
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
        byte[] sequence            = "TACTGACTTTTCTCGCTATTCGTATGCATG".getBytes();
        byte[] expectedOrientation = "CATGCATACGAATAGCGAGAAAAGTCAGTA".getBytes();

        byte[] computedOrientation = SequenceUtils.alphanumericallyLowestOrientation(sequence);

        Assert.assertEquals(expectedOrientation, computedOrientation);
    }

    @Test
    public void computeN50() {
        List<String> sequences = new ArrayList<String>();
        sequences.add("AAGCTTA");
        sequences.add("TTGA");
        sequences.add("AAC");
        sequences.add("TT");
        sequences.add("AA");
        sequences.add("C");
        sequences.add("G");

        int n50 = SequenceUtils.computeN50Value(sequences);

        Assert.assertEquals(4, n50);

        List<CortexKmer> csequences = new ArrayList<CortexKmer>();
        csequences.add(new CortexKmer("AAGCTTA"));
        csequences.add(new CortexKmer("TTGA"));
        csequences.add(new CortexKmer("AAC"));
        csequences.add(new CortexKmer("TT"));
        csequences.add(new CortexKmer("AA"));
        csequences.add(new CortexKmer("C"));
        csequences.add(new CortexKmer("G"));

        int cn50 = SequenceUtils.computeN50Value(csequences);

        Assert.assertEquals(4, cn50);
    }

    /*
    @Test
    public void testEditDistanceGeneration() {
        //byte[] template = "TAC".getBytes();
        byte[] template = SequenceUtils.generateRandomNucleotideSequenceOfLengthN(31);

        Collection<byte[]> seqs1 = SequenceUtils.generateSequencesWithEditDistance1(template);

        System.out.println(new String(template));
        for (byte[] seq : seqs1) {
            System.out.println(new String(seq) + " " + SequenceUtils.editDistance(template, seq));
        }
        System.out.println(seqs1.size());

        Collection<byte[]> seqs2 = SequenceUtils.generateSequencesWithEditDistance2(template);

        System.out.println();
        System.out.println(new String(template));
        for (byte[] seq : seqs2) {
            System.out.println(new String(seq) + " " + SequenceUtils.editDistance(template, seq));
        }
        System.out.println(seqs2.size());
    }
    */
}
