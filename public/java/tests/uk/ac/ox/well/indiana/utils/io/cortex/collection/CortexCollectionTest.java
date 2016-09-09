package uk.ac.ox.well.indiana.utils.io.cortex.collection;

import org.junit.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.File;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

public class CortexCollectionTest {
    @Test
    public void testDynamicGraphMerging() {
        CortexCollection cc = new CortexCollection("testdata/graph.list.txt");

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
        CortexCollection cc = new CortexCollection("testdata/two_samples.sorted.ctx");

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

    @Test
    public void testIteration() {
        CortexCollection cc = new CortexCollection("testdata/graph.list.txt");

        CortexGraph g0 = new CortexGraph(cc.getGraph(0).getCortexFile());
        CortexGraph g1 = new CortexGraph(cc.getGraph(1).getCortexFile());

        Set<String> kmers = new TreeSet<>();
        for (CortexRecord cr : g0) { kmers.add(cr.getKmerAsString()); }
        for (CortexRecord cr : g1) { kmers.add(cr.getKmerAsString()); }

        Iterator<String> it = kmers.iterator();

        for (CortexRecord cr : cc) {
            if (it.hasNext()) {
                String esk = it.next();
                String osk = cr.getKmerAsString();

                Assert.assertEquals(esk, osk);

                CortexRecord c0 = g0.findRecord(esk);
                CortexRecord c1 = g1.findRecord(esk);

                if (c0 != null) {
                    Assert.assertEquals(c0.getCoverage(0), cr.getCoverage(0));
                    Assert.assertEquals(c0.getEdgesAsString(0), cr.getEdgesAsString(0));
                }

                if (c1 != null) {
                    Assert.assertEquals(c1.getCoverage(0), cr.getCoverage(1));
                    Assert.assertEquals(c1.getEdgesAsString(0), cr.getEdgesAsString(1));
                }
            }
        }
    }
}
