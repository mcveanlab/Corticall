package uk.ac.ox.well.indiana.utils.sequence;

import org.testng.Assert;
import org.testng.annotations.Test;

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
    }
}
