package uk.ac.ox.well.cortexjdk.utils.io.graph.links;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import uk.ac.ox.well.cortexjdk.utils.assembler.TempGraphAssembler;
import uk.ac.ox.well.cortexjdk.utils.assembler.TempLinksAssembler;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;

import java.io.File;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;

public class CortexLinksTest {
    @DataProvider(name = "constructLinkData")
    private Object[][] constructLinkData() {
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("test", Collections.singletonList("ACTGATTTCGATGCGATGCGATGCCACGGTGG"));

        Map<String, Collection<String>> reads = new LinkedHashMap<>();
        reads.put("test", Collections.singletonList("TTTCGATGCGATGCGATGCCACG"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 5);
        CortexLinks l = TempLinksAssembler.buildLinks(g, reads, "test");

        return new Object[][] {{l.getFile()}};
    }

    @Test(dataProvider = "constructLinkData")
    public void testLoad(File linksFile) {
        CortexLinksIterable ctp = new CortexLinksIterable(linksFile);

        Assert.assertEquals(4, ctp.getVersion());
        Assert.assertEquals(1, ctp.getNumColors());
        Assert.assertEquals(5, ctp.getKmerSize());
        Assert.assertEquals(21, ctp.getNumKmersInGraph());
        Assert.assertEquals(4, ctp.getNumKmersWithLinks());
        Assert.assertEquals(4, ctp.getNumLinks());

        int numKmersWithLinks = 0;
        int numLinks = 0;
        for (CortexLinksRecord cpr : ctp) {
            numKmersWithLinks++;
            numLinks += cpr.getJunctions().size();
        }

        Assert.assertEquals(4, numKmersWithLinks);
        Assert.assertEquals(ctp.getNumLinks(), numLinks);
    }

    @Test(dataProvider = "constructLinkData")
    public void testMultipleIteration(File linksFile) {
        CortexLinksIterable ctp = new CortexLinksIterable(linksFile);

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
