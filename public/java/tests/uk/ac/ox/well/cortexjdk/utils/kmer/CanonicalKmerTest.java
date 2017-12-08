package uk.ac.ox.well.cortexjdk.utils.kmer;

import org.testng.Assert;
import org.testng.annotations.Test;

public class CanonicalKmerTest {
    @Test
    public void hashcodesAreEqualObjectsAreNot() {
        CanonicalKmer ck1 = new CanonicalKmer("GAACAAAAAAACTTGATAAATGTTTACAAAA");
        CanonicalKmer ck2 = new CanonicalKmer("ACTCTTTTTTAAATGATTATTGCAGATATAT");

        Assert.assertEquals(ck1.hashCode(), ck2.hashCode());
        Assert.assertNotEquals(ck1, ck2);
    }
}
