package uk.ac.ox.well.indiana.utils.sequence;

import org.testng.Assert;
import org.testng.annotations.Test;

public class SequenceUtilsTest {
    @Test
    public void getReverseComplementTest() {
        byte[] sequence   = "TACTGACTTTTCTCGCTATTCGTATGCATG".getBytes();
        byte[] expectedRC = "CATGCATACGAATAGCGAGAAAAGTCAGTA".getBytes();

        byte[] computedRC = SequenceUtils.getReverseComplement(sequence);

        Assert.assertEquals(expectedRC, computedRC);
    }

    @Test
    public void getReverseComplementTestWithNs() {
        byte[] sequence   = "NACTGACTTTTCTCGCTATTCGTATGCATG".getBytes();
        byte[] expectedRC = "CATGCATACGAATAGCGAGAAAAGTCAGTN".getBytes();

        byte[] computedRC = SequenceUtils.getReverseComplement(sequence);

        Assert.assertEquals(expectedRC, computedRC);
    }

    @Test
    public void getReverseComplementWithLowercaseBases() {
        byte[] sequence   = "NACTGACTTTTCTCGCTATTCGTATGCATg".getBytes();
        byte[] expectedRC = "cATGCATACGAATAGCGAGAAAAGTCAGTN".getBytes();

        byte[] computedRC = SequenceUtils.getReverseComplement(sequence);

        Assert.assertEquals(expectedRC, computedRC);
    }

    @Test
    public void getAlphanumericallyLowestOrientation() {
        byte[] sequence            = "TACTGACTTTTCTCGCTATTCGTATGCATG".getBytes();
        byte[] expectedOrientation = "CATGCATACGAATAGCGAGAAAAGTCAGTA".getBytes();

        byte[] computedOrientation = SequenceUtils.getAlphanumericallyLowestOrientation(sequence);

        Assert.assertEquals(expectedOrientation, computedOrientation);
    }
}
