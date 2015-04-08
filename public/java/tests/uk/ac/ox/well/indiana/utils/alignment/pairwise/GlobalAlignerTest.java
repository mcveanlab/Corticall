package uk.ac.ox.well.indiana.utils.alignment.pairwise;

import org.apache.commons.math3.util.Pair;
import org.testng.Assert;
import org.testng.annotations.Test;

public class GlobalAlignerTest {
    @Test
    public void testGlobalAligner() {
        String target = "TTACG";
        String query = "TAG";

        GlobalAligner ga = new GlobalAligner();
        Pair<String, String> a = ga.align(query, target);

        Assert.assertEquals("T-A-G", a.getFirst());
        Assert.assertEquals("TTACG", a.getSecond());
    }
}
