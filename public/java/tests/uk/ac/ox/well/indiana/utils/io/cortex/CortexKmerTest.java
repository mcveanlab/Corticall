package uk.ac.ox.well.indiana.utils.io.cortex;

import org.testng.Assert;
import org.testng.annotations.Test;

public class CortexKmerTest {
    @Test
    public void hashcodesAreEqualObjectsAreNot() {
        CortexKmer ck1 = new CortexKmer("GAACAAAAAAACTTGATAAATGTTTACAAAA");
        CortexKmer ck2 = new CortexKmer("ACTCTTTTTTAAATGATTATTGCAGATATAT");

        Assert.assertEquals(ck1.hashCode(), ck2.hashCode());
        Assert.assertNotEquals(ck1, ck2);
    }
}
