package uk.ac.ox.well.indiana.utils.io.cortex.paths;

import org.testng.Assert;
import org.testng.annotations.Test;

public class CortexPathsTest {
    @Test
    public void testLoadCortexPaths() {
        CortexPaths ctp = new CortexPaths("testdata/PG0085-C.infer.se.ctp");

        Assert.assertEquals(2, ctp.getVersion());
        Assert.assertEquals(1, ctp.getNumColors());
        Assert.assertEquals(47, ctp.getKmerSize());
        Assert.assertEquals(24062238, ctp.getNumKmersInGraph());
        Assert.assertEquals(362383, ctp.getNumKmersWithPaths());
        Assert.assertEquals(2753315, ctp.getNumPaths());
        Assert.assertEquals(7378352, ctp.getPathBytes());

        int numKmersWithPaths = 0;
        int numPaths = 0;
        for (CortexPathsRecord cpr : ctp) {
            numKmersWithPaths++;
            numPaths += cpr.getJunctions().size();
        }

        Assert.assertEquals(ctp.getNumKmersWithPaths(), numKmersWithPaths);
        Assert.assertEquals(ctp.getNumPaths(), numPaths);
    }

    @Test
    public void testMultipleIteration() {
        CortexPaths ctp = new CortexPaths("testdata/PG0085-C.infer.se.ctp");

        CortexPathsRecord cprFirst1 = null, cprLast1 = null;
        for (CortexPathsRecord cpr : ctp) {
            if (cprFirst1 == null) { cprFirst1 = cpr; }
            cprLast1 = cpr;
        }

        CortexPathsRecord cprFirst2 = null, cprLast2 = null;
        for (CortexPathsRecord cpr : ctp) {
            if (cprFirst2 == null) { cprFirst2 = cpr; }
            cprLast2 = cpr;
        }

        Assert.assertEquals(cprFirst1, cprFirst2);
        Assert.assertEquals(cprLast1, cprLast2);
        Assert.assertNotEquals(cprFirst1, cprLast1);
    }
}
