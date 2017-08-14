package uk.ac.ox.well.cortexjdk.utils.io.cortex.links;

import org.testng.Assert;
import org.testng.annotations.Test;

public class CortexLinksTest {
    @Test
    public void testLoadFormatVersion2() {
        CortexLinksIterable ctp = new CortexLinksIterable("testdata/PG0085-C.infer.se.ctp");

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
    public void testLoadFormatVersion3() {
        CortexLinksIterable ctp = new CortexLinksIterable("testdata/PG0063-C.ERR019060.infer.pe.k51.v3.ctp");

        Assert.assertEquals(3, ctp.getVersion());
        Assert.assertEquals(1, ctp.getNumColors());
        Assert.assertEquals(51, ctp.getKmerSize());
        Assert.assertEquals(21826964, ctp.getNumKmersInGraph());
        Assert.assertEquals(71313, ctp.getNumKmersWithLinks());
        Assert.assertEquals(1061391, ctp.getNumLinks());
        Assert.assertEquals(2511987, ctp.getLinkBytes());

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
        CortexLinksIterable ctp = new CortexLinksIterable("testdata/PG0085-C.infer.se.ctp");

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

    /*
    @Test
    public void testParseRecordWithExtendedInfo() {
        CortexGraphLinks ctp = new CortexGraphLinks("testdata/PG0051-C.ERR019061.chr1.se.ctp");

        for (CortexLinksRecord clr : ctp) {
            for (CortexJunctionsRecord cjr : clr.getJunctions()) {
                String seq = cjr.getSeq();

                Assert.assertNotNull(seq);
            }
        }

        CortexGraphLinks ctp2 = new CortexGraphLinks("testdata/PG0085-C.infer.se.ctp");

        for (CortexLinksRecord clr : ctp2) {
            for (CortexJunctionsRecord cjr : clr.getJunctions()) {
                String seq = cjr.getSeq();

                Assert.assertNull(seq);
            }
        }
    }
    */
}
