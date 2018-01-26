package uk.ac.ox.well.cortexjdk.utils.io.graph.collection;

import org.junit.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import uk.ac.ox.well.cortexjdk.utils.assembler.TempGraphAssembler;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexCollection;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.util.*;

public class CortexCollectionTest {
    @DataProvider(name="joinedTempGraph")
    public static Object[][] createJoinedTempGraph() {
        Map<String, Collection<String>> h1 = new LinkedHashMap<>();
        h1.put("fred", Collections.singletonList("ACCGTATGTA"));
        h1.put("wilma", Collections.singletonList("ACCGCATGTA"));
        h1.put("barney", Collections.singletonList("ACCGTATATA"));

        Map<String, Collection<String>> h2 = new LinkedHashMap<>();
        h2.put("pebbles", Collections.singletonList("AACGTATGTA"));
        h2.put("bam-bam", Collections.singletonList("ACTGCATGTA"));

        CortexGraph g1 = TempGraphAssembler.buildGraph(h1, 3);
        CortexGraph g2 = TempGraphAssembler.buildGraph(h2, 3);

        CortexCollection cc = new CortexCollection(g1, g2);

        return new Object[][] {{cc, g1, g2}};
    }

    @Test(dataProvider = "joinedTempGraph")
    public void testDynamicGraphMerging(CortexCollection cc, CortexGraph g0, CortexGraph g1) {
        for (CortexRecord c0 : g0) {
            CortexRecord c1 = g1.findRecord(c0.getKmerAsString());
            CortexRecord cm = cc.findRecord(c0.getKmerAsString());

            Assert.assertEquals(c0.getKmerAsString(), cm.getKmerAsString());

            Assert.assertEquals(c0.getCoverage(0), cm.getCoverage(0));
            Assert.assertEquals(c0.getCoverage(1), cm.getCoverage(1));
            Assert.assertEquals(c0.getCoverage(2), cm.getCoverage(2));

            Assert.assertEquals(c0.getEdgesAsString(0), cm.getEdgesAsString(0));
            Assert.assertEquals(c0.getEdgesAsString(1), cm.getEdgesAsString(1));
            Assert.assertEquals(c0.getEdgesAsString(2), cm.getEdgesAsString(2));

            if (c1 == null) {
                Assert.assertEquals(0, cm.getCoverage(3));
                Assert.assertEquals(0, cm.getCoverage(4));
                Assert.assertEquals("........", cm.getEdgesAsString(3));
                Assert.assertEquals("........", cm.getEdgesAsString(4));
            } else {
                Assert.assertEquals(c1.getCoverage(0), cm.getCoverage(3));
                Assert.assertEquals(c1.getCoverage(1), cm.getCoverage(4));

                Assert.assertEquals(c1.getEdgesAsString(0), cm.getEdgesAsString(3));
                Assert.assertEquals(c1.getEdgesAsString(1), cm.getEdgesAsString(4));
            }
        }
    }

    @Test(dataProvider = "joinedTempGraph")
    public void testSingleFileCollection(CortexCollection ccUnused, CortexGraph g0, CortexGraph g1Unused) {
        CortexCollection cc = new CortexCollection(g0);

        for (CortexRecord c0 : g0) {
            CortexRecord cm = cc.findRecord(c0.getKmerAsString());

            Assert.assertEquals(c0.getKmerAsString(), cm.getKmerAsString());

            Assert.assertEquals(c0.getCoverage(0), cm.getCoverage(0));
            Assert.assertEquals(c0.getCoverage(1), cm.getCoverage(1));

            Assert.assertEquals(c0.getEdgesAsString(0), cm.getEdgesAsString(0));
            Assert.assertEquals(c0.getEdgesAsString(1), cm.getEdgesAsString(1));
        }
    }

    @Test(dataProvider = "joinedTempGraph")
    public void testIteration(CortexCollection cc, CortexGraph g0, CortexGraph g1) {
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

    @Test(dataProvider = "joinedTempGraph")
    public void testNumColors(CortexCollection cc, CortexGraph g1, CortexGraph g2) {
        Assert.assertEquals(g1.getNumColors() + g2.getNumColors(), cc.getNumColors());
    }

    @Test(dataProvider = "joinedTempGraph")
    public void testColorNames(CortexCollection cc, CortexGraph g1, CortexGraph g2) {
        Assert.assertEquals(g1.getSampleName(0), cc.getSampleName(0));
        Assert.assertEquals(g1.getSampleName(1), cc.getSampleName(1));
        Assert.assertEquals(g1.getSampleName(2), cc.getSampleName(2));
        Assert.assertEquals(g2.getSampleName(0), cc.getSampleName(3));
        Assert.assertEquals(g2.getSampleName(1), cc.getSampleName(4));
    }

    @Test(dataProvider = "joinedTempGraph")
    public void testColors(CortexCollection cc, CortexGraph g1, CortexGraph g2) {
        Assert.assertEquals(g1.getColor(0), cc.getColor(0));
        Assert.assertEquals(g1.getColor(1), cc.getColor(1));
        Assert.assertEquals(g1.getColor(2), cc.getColor(2));
        Assert.assertEquals(g2.getColor(0), cc.getColor(3));
        Assert.assertEquals(g2.getColor(1), cc.getColor(4));
    }
}
