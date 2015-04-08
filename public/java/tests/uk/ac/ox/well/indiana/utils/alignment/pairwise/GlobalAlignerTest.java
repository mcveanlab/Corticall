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

        Assert.assertEquals(a.getFirst(), "T-A-G");
        Assert.assertEquals(a.getSecond(), "TTACG");
    }
}
