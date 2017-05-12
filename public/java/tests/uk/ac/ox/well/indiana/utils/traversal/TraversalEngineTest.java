package uk.ac.ox.well.indiana.utils.traversal;

import org.jetbrains.annotations.NotNull;
import org.jgrapht.DirectedGraph;
import org.jgrapht.traverse.TopologicalOrderIterator;
import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.*;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.stoppingconditions.LongWalkStopper;

import java.io.File;
import java.io.IOException;
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

    private CortexGraph buildGraph(Map<String, String> haplotypes, int kmerSize) {
        File tempFile = null;
        try {
            tempFile = File.createTempFile("tempgraph", ".ctx");
        } catch (IOException e) {
            throw new IndianaException("Could not get a temp file for graph creation");
        }
        tempFile.deleteOnExit();

        CortexGraphWriter cgw = new CortexGraphWriter(tempFile);
        cgw.setHeader(constructCortexHeader(haplotypes, kmerSize));

        Map<CortexKmer, CortexRecord> crs = new TreeMap<>();

        int c = 0;
        for (String sampleName : haplotypes.keySet()) {
            String sequence = haplotypes.get(sampleName);

            for (int i = 0; i <= sequence.length() - kmerSize; i++) {
                String sk = sequence.substring(i, i + kmerSize);

                String prevBase = i == 0 ? null : sequence.substring(i - 1, i);
                String nextBase = i == sequence.length() - kmerSize ? null : sequence.substring(i + kmerSize, i + kmerSize + 1);

                updateRecord(crs, haplotypes.size(), c, sk, prevBase, nextBase);
            }

            c++;
        }

        for (CortexKmer cr : crs.keySet()) {
            cgw.addRecord(crs.get(cr));
        }

        cgw.close();

        return new CortexGraph(tempFile);
    }

    private void updateRecord(Map<CortexKmer, CortexRecord> crs, int numColors, int color, String sk, String prevBase, String nextBase) {
        CortexKmer ck = new CortexKmer(sk);

        List<Integer> coverageList = new ArrayList<>();
        List<Set<String>> inEdgesList = new ArrayList<>();
        List<Set<String>> outEdgesList = new ArrayList<>();

        CortexRecord oldRc = crs.containsKey(ck) ? crs.get(ck) : null;

        for (int c = 0; c < numColors; c++) {
            int cov = 0;
            Set<String> inEdges = new HashSet<>();
            Set<String> outEdges = new HashSet<>();

            if (oldRc != null) {
                cov += oldRc.getCoverage(c);

                inEdges.addAll(!ck.isFlipped() ? oldRc.getInEdgesAsStrings(c, false) : oldRc.getOutEdgesAsStrings(c, true));
                outEdges.addAll(!ck.isFlipped() ? oldRc.getOutEdgesAsStrings(c, false) : oldRc.getInEdgesAsStrings(c, true));
            }

            if (c == color) {
                cov++;

                if (!ck.isFlipped()) {
                    if (prevBase != null) { inEdges.add(prevBase); }
                    if (nextBase != null) { outEdges.add(nextBase); }
                } else {
                    if (nextBase != null) { inEdges.add(SequenceUtils.reverseComplement(nextBase)); }
                    if (prevBase != null) { outEdges.add(SequenceUtils.reverseComplement(prevBase)); }
                }
            }

            coverageList.add(cov);
            inEdgesList.add(inEdges);
            outEdgesList.add(outEdges);
        }

        crs.put(ck, new CortexRecord(SequenceUtils.alphanumericallyLowestOrientation(sk), coverageList, inEdgesList, outEdgesList));
    }

    @NotNull
    private CortexHeader constructCortexHeader(Map<String, String> haplotypes, int kmerSize) {
        CortexHeader ch = new CortexHeader();
        ch.setVersion(6);
        ch.setNumColors(haplotypes.size());
        ch.setKmerSize(kmerSize);
        ch.setKmerBits(CortexRecord.getKmerBits(kmerSize));

        for (String sampleName : haplotypes.keySet()) {
            CortexColor col = new CortexColor();

            col.setSampleName(sampleName);
            col.setCleanedAgainstGraph(false);
            col.setCleanedAgainstGraphName("");
            col.setErrorRate(0.0);
            col.setLowCovgKmersRemoved(false);
            col.setLowCovgSupernodesRemoved(false);
            col.setLowCovKmerThreshold(0);
            col.setTipClippingApplied(false);
            col.setMeanReadLength(0);
            col.setTotalSequence(0);

            ch.addColor(col);
        }

        return ch;
    }

    @Test
    public void testArbitraryGraphConstruction() {
        Map<String, String> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom", "AATA");
        haplotypes.put("dad", "AATG");

        CortexGraph g = buildGraph(haplotypes, 3);

        CortexGraphExpectation cge = new CortexGraphExpectation(g);

        cge.hasNRecords(3);
        cge.hasRecord("AAT 1 1 ....A... ......G.");
        cge.hasRecord("ATA 1 0 a....... ........");
        cge.hasRecord("ATG 0 1 ........ a.......");
    }

    @Test
    public void testTraversalEngine() {
        CortexGraph g = new CortexGraph("3D7xHB3.k47.clean.infer.ctx");

        int childColor = g.getColorForSampleName("PG0051-C.ERR019061.md");
        Set<Integer> parentColors = new HashSet<>(g.getColorsForSampleNames(Arrays.asList("PG0051-C.ERR019061.md", "PG0052-C.ERR019054.md")));
        Set<Integer> refColors = new HashSet<>(g.getColorsForSampleNames(Arrays.asList("3D7", "HB3")));

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(childColor)
                .joiningColors(parentColors)
                .recruitmentColors(refColors)
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.AND)
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
                .connectUnusedNeighbors(false)
                .stopper(new LongWalkStopper())
                .graph(g)
                .make();

        CortexKmer ck = new CortexKmer("ATGGAAATTGTATAAATAAAGATAATGACAACACATGTAAAAATTCA");
        String sk = ck.getKmerAsString();

        DirectedGraph<CortexVertex, CortexEdge> walkNew = e.dfs(sk);

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

        System.out.println(sb);
        System.out.println();
    }
}
