package uk.ac.ox.well.indiana.utils.traversal;

import org.jgrapht.DirectedGraph;
import org.jgrapht.traverse.TopologicalOrderIterator;
import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.assembler.TempGraphAssembler;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.*;
import uk.ac.ox.well.indiana.utils.stoppingconditions.ContigStopper;

import java.util.*;

/**
 * Created by kiran on 10/05/2017.
 */
public class TraversalEngineTest {
    private class CortexGraphExpectation {
        private CortexGraph eval;

        public CortexGraphExpectation(CortexGraph eval) {
            this.eval = eval;
        }

        public void hasNRecords(int n) {
            Assert.assertEquals(eval.getNumRecords(), n);
        }

        public void hasRecord(String truth) {
            boolean hasRecord = false;
            for (CortexRecord cr : eval) {
                if (cr.toString().equals(truth)) {
                    hasRecord = true;
                    break;
                }
            }

            Assert.assertTrue(hasRecord, "Did not find record '" + truth + "'");
        }
    }

    @Test
    public void testArbitraryGraphConstruction() {
        Map<String, String> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom", "AATA");
        haplotypes.put("dad", "AATG");

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 3);

        CortexGraphExpectation cge = new CortexGraphExpectation(g);

        cge.hasNRecords(3);
        cge.hasRecord("AAT 1 1 ....A... ......G.");
        cge.hasRecord("ATA 1 0 a....... ........");
        cge.hasRecord("ATG 0 1 ........ a.......");
    }

    @Test
    public void testSlightlyLargerArbitraryGraphConstruction() {
        Map<String, String> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom", "AGTTCTGATCTGGGCTATATGCT");
        haplotypes.put("dad", "AGTTCTGATCTGGGCTATATGCT");
        haplotypes.put("kid", "AGTTCTGATCTGGGCTATATGCT");

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 5);

        CortexGraphExpectation cge = new CortexGraphExpectation(g);

        cge.hasNRecords(19);
        cge.hasRecord("AGAAC 1 1 1 .c.....T .c.....T .c.....T");
        cge.hasRecord("AGATC 1 1 1 .c..A... .c..A... .c..A...");
        cge.hasRecord("AGCAT 1 1 1 ....A... ....A... ....A...");
        cge.hasRecord("AGCCC 1 1 1 ...tA... ...tA... ...tA...");
        cge.hasRecord("AGTTC 1 1 1 .......T .......T .......T");
        cge.hasRecord("ATAGC 1 1 1 ...t.C.. ...t.C.. ...t.C..");
        cge.hasRecord("ATATA 1 1 1 .c....G. .c....G. .c....G.");
        cge.hasRecord("ATATG 1 1 1 ...t.C.. ...t.C.. ...t.C..");
        cge.hasRecord("ATCAG 1 1 1 ..g.A... ..g.A... ..g.A...");
        cge.hasRecord("ATCTG 1 1 1 ..g...G. ..g...G. ..g...G.");
        cge.hasRecord("CAGAA 1 1 1 ...t.C.. ...t.C.. ...t.C..");
        cge.hasRecord("CCAGA 1 1 1 .c.....T .c.....T .c.....T");
        cge.hasRecord("CCCAG 1 1 1 ..g.A... ..g.A... ..g.A...");
        cge.hasRecord("CTATA 1 1 1 ..g....T ..g....T ..g....T");
        cge.hasRecord("GATCA 1 1 1 a.....G. a.....G. a.....G.");
        cge.hasRecord("GCATA 1 1 1 a......T a......T a......T");
        cge.hasRecord("GCCCA 1 1 1 a.....G. a.....G. a.....G.");
        cge.hasRecord("GGCTA 1 1 1 ..g....T ..g....T ..g....T");
        cge.hasRecord("TCAGA 1 1 1 a...A... a...A... a...A...");
    }

    @Test
    public void testShortContigReconstruction() {
        Map<String, String> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom", "AGTTCTGATCTGGGCTATATGCT");
        haplotypes.put("dad", "AGTTCGAATCTGGGCTATATGCT");
        haplotypes.put("kid", "AGTTCTGATCTGGGCTATGGCTA");

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 5);

        TraversalEngine e = new TraversalEngineFactory()
                .joiningColors(g.getColorsForSampleNames(Arrays.asList("mom", "dad")))
                .recruitmentColors(g.getColorsForSampleNames(Arrays.asList("mom", "dad")))
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.AND)
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
                .connectUnusedNeighbors(false)
                .stopper(new ContigStopper())
                .graph(g)
                .make();

        Map<String, String> expectations = new LinkedHashMap<>();
        expectations.put("mom", "AGTTCTGATCTGGGCTATATGCT");
        expectations.put("dad",   "TTCGAATCTGGGCTATATGCT");
        expectations.put("kid", "AGTTCTGATCTGGGCTATGGCT");

        for (int c = 0; c < 3; c++) {
            e.getConfiguration().setTraversalColor(c);

            DirectedGraph<CortexVertex, CortexEdge> walkNew = e.dfs("CTGGG");

            TopologicalOrderIterator<CortexVertex, CortexEdge> toi = new TopologicalOrderIterator<>(walkNew);

            StringBuilder sb = new StringBuilder();
            while (toi.hasNext()) {
                CortexVertex cv = toi.next();
                if (sb.length() == 0) {
                    sb.append(cv.getSk());
                } else {
                    sb.append(cv.getSk().substring(cv.getSk().length() - 1));
                }
            }

            Assert.assertEquals(sb.toString(), expectations.get(g.getSampleName(c)));
        }
    }

    @Test
    public void testTraversalEngine() {
        Map<String, String> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom", "AGTTCTGATCTGGGCTATATGCT");
        haplotypes.put("dad", "AGTTCGAATCTGGGCTATATGCT");
        haplotypes.put("kid", "AGTTCTGATCTGGGCTATGGCTA");

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 5);

        Set<Integer> parentColors = new HashSet<>(g.getColorsForSampleNames(Arrays.asList("mom", "dad")));
        Set<Integer> refColors = new HashSet<>(g.getColorsForSampleNames(Arrays.asList("mom", "dad")));

        TraversalEngine e = new TraversalEngineFactory()
                .joiningColors(parentColors)
                .recruitmentColors(refColors)
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.AND)
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
                .connectUnusedNeighbors(false)
                .stopper(new ContigStopper())
                .graph(g)
                .make();

        Map<String, String> expectations = new LinkedHashMap<>();
        expectations.put("mom", "AGTTCTGATCTGGGCTATATGCT");
        expectations.put("dad",   "TTCGAATCTGGGCTATATGCT");
        expectations.put("kid", "AGTTCTGATCTGGGCTATGGCT" );

        for (int c = 0; c < 3; c++) {
            e.getConfiguration().setTraversalColor(c);

            DirectedGraph<CortexVertex, CortexEdge> walkNew = e.dfs("CTGGG");

            TopologicalOrderIterator<CortexVertex, CortexEdge> toi = new TopologicalOrderIterator<>(walkNew);

            StringBuilder sb = new StringBuilder();

            while (toi.hasNext()) {
                CortexVertex cv = toi.next();
                if (sb.length() == 0) {
                    sb.append(cv.getSk());
                } else {
                    sb.append(cv.getSk().substring(cv.getSk().length() - 1));
                }
            }

            Assert.assertEquals(sb.toString(), expectations.get(g.getSampleName(c)));
        }
    }
}
