package uk.ac.ox.well.indiana.utils.io.cortex.collection;

import org.junit.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.File;

public class CortexCollectionTest {
    @Test
    public void testDynamicGraphMerging() {
        CortexCollection cc = new CortexCollection(new File("testdata/graph.list.txt"));

        CortexGraph g0 = new CortexGraph(cc.getGraph(0).getCortexFile());
        CortexGraph g1 = new CortexGraph(cc.getGraph(1).getCortexFile());

        for (CortexRecord c0 : g0) {
            CortexRecord c1 = g1.findRecord(c0.getKmerAsString());
            CortexRecord cm = cc.findRecord(c0.getKmerAsString());

            Assert.assertEquals(c0.getKmerAsString(), cm.getKmerAsString());

            Assert.assertEquals(c0.getCoverage(0), cm.getCoverage(0));
            Assert.assertEquals(c0.getEdgesAsString(0), cm.getEdgesAsString(0));

            if (c1 == null) {
                Assert.assertEquals(0, cm.getCoverage(1));
                Assert.assertEquals("........", cm.getEdgesAsString(1));
            } else {
                Assert.assertEquals(c1.getCoverage(0), cm.getCoverage(1));
                Assert.assertEquals(c1.getEdgesAsString(0), cm.getEdgesAsString(1));
            }
        }
    }

    @Test
    public void testSingleFileCollection() {
        CortexCollection cc = new CortexCollection(new File("testdata/two_samples.sorted.ctx"));

        Assert.assertEquals(2, cc.getNumColors());

        CortexGraph g0 = new CortexGraph(cc.getGraph(0).getCortexFile());

        for (CortexRecord c0 : g0) {
            CortexRecord cm = cc.findRecord(c0.getKmerAsString());

            Assert.assertEquals(c0.getKmerAsString(), cm.getKmerAsString());

            Assert.assertEquals(c0.getCoverage(0), cm.getCoverage(0));
            Assert.assertEquals(c0.getCoverage(1), cm.getCoverage(1));

            Assert.assertEquals(c0.getEdgesAsString(0), cm.getEdgesAsString(0));
            Assert.assertEquals(c0.getEdgesAsString(1), cm.getEdgesAsString(1));
        }
    }
}
