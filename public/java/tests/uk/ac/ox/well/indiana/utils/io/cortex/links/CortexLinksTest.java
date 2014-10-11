package uk.ac.ox.well.indiana.utils.io.cortex.links;

import org.testng.Assert;
import org.testng.annotations.Test;

public class CortexLinksTest {
    @Test
    public void testLoadCortexLinks() {
        CortexLinks ctp = new CortexLinks("testdata/PG0085-C.infer.se.ctp");

        Assert.assertEquals(2, ctp.getVersion());
        Assert.assertEquals(1, ctp.getNumColors());
        Assert.assertEquals(47, ctp.getKmerSize());
        Assert.assertEquals(24062238, ctp.getNumKmersInGraph());
        Assert.assertEquals(362383, ctp.getNumKmersWithLinks());
        Assert.assertEquals(2753315, ctp.getNumLinks());
        Assert.assertEquals(7378352, ctp.getLinkBytes());

        int numKmersWithLinks = 0;
        int numLinks = 0;
        for (CortexLinksRecord cpr : ctp) {
            numKmersWithLinks++;
            numLinks += cpr.getJunctions().size();
        }

        Assert.assertEquals(ctp.getNumKmersWithLinks(), numKmersWithLinks);
        Assert.assertEquals(ctp.getNumLinks(), numLinks);
    }

    @Test
    public void testMultipleIteration() {
        CortexLinks ctp = new CortexLinks("testdata/PG0085-C.infer.se.ctp");

        CortexLinksRecord cprFirst1 = null, cprLast1 = null;
        for (CortexLinksRecord cpr : ctp) {
            if (cprFirst1 == null) { cprFirst1 = cpr; }
            cprLast1 = cpr;
        }

        CortexLinksRecord cprFirst2 = null, cprLast2 = null;
        for (CortexLinksRecord cpr : ctp) {
            if (cprFirst2 == null) { cprFirst2 = cpr; }
            cprLast2 = cpr;
        }

        Assert.assertEquals(cprFirst1, cprFirst2);
        Assert.assertEquals(cprLast1, cprLast2);
        Assert.assertNotEquals(cprFirst1, cprLast1);
    }
}
