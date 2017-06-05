package uk.ac.ox.well.indiana.utils.traversal;

import htsjdk.samtools.reference.FastaSequenceFile;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.DirectedGraph;
import org.jgrapht.GraphPath;
import org.jgrapht.Graphs;
import org.jgrapht.alg.DijkstraShortestPath;
import org.jgrapht.alg.shortestpath.KShortestPaths;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.assembler.TempGraphAssembler;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.*;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksMap;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.stoppingconditions.ContigStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.CycleCollapsingContigStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.DestinationStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.ExplorationStopper;
import uk.ac.ox.well.indiana.utils.visualizer.GraphVisualizer;

import java.io.File;
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
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom", Collections.singletonList("AATA"));
        haplotypes.put("dad", Collections.singletonList("AATG"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 3);

        CortexGraphExpectation cge = new CortexGraphExpectation(g);

        cge.hasNRecords(3);
        cge.hasRecord("AAT 1 1 ....A... ......G.");
        cge.hasRecord("ATA 1 0 a....... ........");
        cge.hasRecord("ATG 0 1 ........ a.......");
    }

    @Test
    public void testSlightlyLargerArbitraryGraphConstruction() {
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom", Collections.singletonList("AGTTCTGATCTGGGCTATATGCT"));
        haplotypes.put("dad", Collections.singletonList("AGTTCTGATCTGGGCTATATGCT"));
        haplotypes.put("kid", Collections.singletonList("AGTTCTGATCTGGGCTATATGCT"));

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
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom", Collections.singletonList("AGTTCTGATCTGGGCTATATGCT"));
        haplotypes.put("dad", Collections.singletonList("AGTTCGAATCTGGGCTATATGCT"));
        haplotypes.put("kid", Collections.singletonList("AGTTCTGATCTGGGCTATGGCTA"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 5);

        TraversalEngine e = new TraversalEngineFactory()
                .joiningColors(g.getColorsForSampleNames(Arrays.asList("mom", "dad")))
                .recruitmentColors(g.getColorsForSampleNames(Arrays.asList("mom", "dad")))
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.AND)
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
                .connectAllNeighbors(false)
                .stopper(ContigStopper.class)
                .graph(g)
                .make();

        Map<String, String> expectations = new LinkedHashMap<>();
        expectations.put("mom", "AGTTCTGATCTGGGCTATATGCT");
        expectations.put("dad",   "TTCGAATCTGGGCTATATGCT");
        expectations.put("kid", "AGTTCTGATCTGGGCTATGGCTA");

        for (int c = 0; c < 3; c++) {
            e.getConfiguration().setTraversalColor(c);

            DirectedWeightedPseudograph<CortexVertex, CortexEdge> walk = e.dfs("CTGGG");

            String contig = e.getContig(walk, "CTGGG", c);

            Assert.assertEquals(contig, expectations.get(g.getSampleName(c)));
        }
    }

    @Test
    public void testRecruitment() {
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom", Collections.singletonList("AGTTCTGATCTGGGCTATATGCT"));
        haplotypes.put("dad", Collections.singletonList("AGTTCTGATCTGGGCTATATGCT"));
        haplotypes.put("kid",             Arrays.asList("AGTTCTG",      "ATGGCTA"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 5);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(g.getColorForSampleName("kid"))
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.AND)
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
                .connectAllNeighbors(false)
                .stopper(ContigStopper.class)
                .graph(g)
                .make();

        Map<Boolean, String> expectations = new HashMap<>();
        expectations.put(true,  "AGTTCTGATCTGGGCTATATGCT");
        expectations.put(false, "AGTTCTG");

        for (boolean useRecruitment : Arrays.asList(true, false)) {
            if (useRecruitment) {
                e.getConfiguration().setRecruitmentColors(g.getColorsForSampleNames(Arrays.asList("mom", "dad")));
            } else {
                e.getConfiguration().setRecruitmentColors();
            }

            DirectedWeightedPseudograph<CortexVertex, CortexEdge> walkNew = e.dfs("GTTCT");

            String contig = e.getContig(walkNew, "GTTCT", g.getColorForSampleName("kid"));

            Assert.assertEquals(contig, expectations.get(useRecruitment));
        }
    }

    @Test
    public void testCyclesAreNotAssembled() {
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("test", Collections.singletonList("GGATCAGTCCAGTCCAGTCCAGTCCCCCCT"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 5);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(g.getColorForSampleName("test"))
                .stopper(CycleCollapsingContigStopper.class)
                .graph(g)
                .make();

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> walk = e.dfs("GATCA");
        String contig = e.getContig(walk, "GATCA", g.getColorForSampleName("test"));

        //Assert.assertEquals("GGATCAGTCCAGTCCCCCCT", contig);
        Assert.assertEquals("GGATCAGTCC", contig);
    }

    @Test
    public void testVarGeneReconstruction() {
        FastaSequenceFile f = new FastaSequenceFile(new File("testdata/var.fasta"), true);
        CortexGraph g = new CortexGraph("testdata/var.ctx");
        CortexLinksMap l = new CortexLinksMap("testdata/var.ctp.gz");

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(g.getColorForSampleName("var"))
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.AND)
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
                .connectAllNeighbors(false)
                .stopper(CycleCollapsingContigStopper.class)
                .graph(g)
                .links(l)
                .make();

        String varSequence = f.nextSequence().getBaseString().toUpperCase();
        String[] kmers = { "ATGTGCCCTTGAATATGAATATTATAAGCATACTAATGGCGGTGGTA", "GTACCAAATGATTATAAAAGTGGTGATATTCCATTGAATACACAACC" };

        for (String kmer : kmers) {
            DirectedWeightedPseudograph<CortexVertex, CortexEdge> walk = e.dfs(kmer);

            String contig = e.getContig(walk, kmer, g.getColorForSampleName("var"));

            Assert.assertEquals(varSequence, contig);
        }
    }

    @Test
    public void testCallSNV() {
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom", Collections.singletonList("AGTTCGAATCTGGGCTATATGCT"));
        haplotypes.put("dad", Collections.singletonList("AGTTCGAATCTGGGCTATATGCT"));
        haplotypes.put("kid", Collections.singletonList("AGTTCGAATCTGCGCTATATGCT"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 7);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(g.getColorForSampleName("kid"))
                .displayColors(g.getColorsForSampleNames(Arrays.asList("mom", "dad", "kid")))
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.AND)
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
                .connectAllNeighbors(false)
                .stopper(ExplorationStopper.class)
                .graph(g)
                .make();

        String seed = "TCTGCGC";

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> walk = e.dfs(seed);

        Set<CortexVertex> iss = new HashSet<>();
        Set<CortexVertex> oss = new HashSet<>();

        for (CortexVertex v : walk.vertexSet()) {
            if (TraversalEngine.inDegree(walk, v) > 1) {
                iss.add(v);
            }

            if (TraversalEngine.outDegree(walk, v) > 1) {
                oss.add(v);
            }
        }

        e.getConfiguration().setStoppingRule(DestinationStopper.class);
        e.getConfiguration().setGraphCombinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.OR);

        for (CortexVertex is : iss) {
            for (CortexVertex os : oss) {
                DirectedGraph<CortexVertex, CortexEdge> p = new DefaultDirectedGraph<>(CortexEdge.class);
                p.addVertex(os);

                e.getConfiguration().setPreviousTraversal(p);

                for (int c : g.getColorsForSampleNames(Arrays.asList("mom", "dad"))) {
                    e.getConfiguration().setTraversalColor(c);

                    DirectedWeightedPseudograph<CortexVertex, CortexEdge> w = e.dfs(is.getSk());

                    if (w != null) {
                        Graphs.addGraph(walk, w);
                    }
                }
            }
        }

        PathFinder dspKid = new PathFinder(walk, g.getColorForSampleName("kid"));
        PathFinder dspMom = new PathFinder(walk, g.getColorForSampleName("mom"));
        PathFinder dspDad = new PathFinder(walk, g.getColorForSampleName("dad"));

        Set<Pair<String, String>> allelesVsMom = new HashSet<>();
        Set<Pair<String, String>> allelesVsDad = new HashSet<>();

        for (CortexVertex os : oss) {
            for (CortexVertex is : iss) {
                GraphPath<CortexVertex, CortexEdge> pk = dspKid.getPathFinder(os, is, new CortexKmer(seed), true);
                GraphPath<CortexVertex, CortexEdge> pm = dspMom.getPathFinder(os, is);
                GraphPath<CortexVertex, CortexEdge> pd = dspDad.getPathFinder(os, is);

                Pair<String, String> akm = pathsToAlleles(pk, pm);
                Pair<String, String> akd = pathsToAlleles(pk, pd);

                allelesVsMom.add(akm);
                allelesVsDad.add(akd);
            }
        }

        Pair<String, String> alleleVsMom = allelesVsMom.iterator().next();
        Pair<String, String> alleleVsDad = allelesVsDad.iterator().next();

        Assert.assertEquals(alleleVsMom.getFirst(), "C");
        Assert.assertEquals(alleleVsMom.getSecond(), "G");
        Assert.assertEquals(alleleVsDad.getFirst(), "C");
        Assert.assertEquals(alleleVsDad.getSecond(), "G");
    }

    @Test
    public void testCallInsertion() {
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom", Collections.singletonList("AGTTCGAATCTGGGCTATATGCT"));
        haplotypes.put("dad", Collections.singletonList("AGTTCGAATCTGGGCTATATGCT"));
        haplotypes.put("kid", Collections.singletonList("AGTTCGAATCTGGCATGCTGCTATATGCT"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 7);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(g.getColorForSampleName("kid"))
                .displayColors(g.getColorsForSampleNames(Arrays.asList("mom", "dad", "kid")))
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.AND)
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
                .connectAllNeighbors(false)
                .stopper(ExplorationStopper.class)
                .graph(g)
                .make();

        String seed = "GCATGCT";

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> walk = e.dfs(seed);

        Set<CortexVertex> iss = new HashSet<>();
        Set<CortexVertex> oss = new HashSet<>();

        for (CortexVertex v : walk.vertexSet()) {
            if (TraversalEngine.inDegree(walk, v) > 1) {
                iss.add(v);
            }

            if (TraversalEngine.outDegree(walk, v) > 1) {
                oss.add(v);
            }
        }

        e.getConfiguration().setStoppingRule(DestinationStopper.class);
        e.getConfiguration().setGraphCombinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.OR);

        for (CortexVertex is : iss) {
            for (CortexVertex os : oss) {
                DirectedGraph<CortexVertex, CortexEdge> p = new DefaultDirectedGraph<>(CortexEdge.class);
                p.addVertex(os);

                e.getConfiguration().setPreviousTraversal(p);

                for (int c : g.getColorsForSampleNames(Arrays.asList("mom", "dad"))) {
                    e.getConfiguration().setTraversalColor(c);

                    DirectedWeightedPseudograph<CortexVertex, CortexEdge> w = e.dfs(is.getSk());

                    if (w != null) {
                        Graphs.addGraph(walk, w);
                    }
                }
            }
        }

        PathFinder dspKid = new PathFinder(walk, g.getColorForSampleName("kid"));
        PathFinder dspMom = new PathFinder(walk, g.getColorForSampleName("mom"));
        PathFinder dspDad = new PathFinder(walk, g.getColorForSampleName("dad"));

        Set<Pair<String, String>> allelesVsMom = new HashSet<>();
        Set<Pair<String, String>> allelesVsDad = new HashSet<>();

        for (CortexVertex os : oss) {
            for (CortexVertex is : iss) {
                GraphPath<CortexVertex, CortexEdge> pk = dspKid.getPathFinder(os, is, new CortexKmer(seed), true);
                GraphPath<CortexVertex, CortexEdge> pm = dspMom.getPathFinder(os, is);
                GraphPath<CortexVertex, CortexEdge> pd = dspDad.getPathFinder(os, is);

                Pair<String, String> akm = pathsToAlleles(pk, pm);
                Pair<String, String> akd = pathsToAlleles(pk, pd);

                allelesVsMom.add(akm);
                allelesVsDad.add(akd);
            }
        }

        Pair<String, String> alleleVsMom = allelesVsMom.iterator().next();
        Pair<String, String> alleleVsDad = allelesVsDad.iterator().next();

        Assert.assertEquals(alleleVsMom.getFirst(), "CATGCT");
        Assert.assertEquals(alleleVsMom.getSecond(), "");
        Assert.assertEquals(alleleVsDad.getFirst(), "CATGCT");
        Assert.assertEquals(alleleVsDad.getSecond(), "");
    }

    @Test
    public void testCallDeletion() {
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom", Collections.singletonList("AGTTCGAATCTGGGCTATATGCT"));
        haplotypes.put("dad", Collections.singletonList("AGTTCGAATCTGGGCTATATGCT"));
        haplotypes.put("kid", Collections.singletonList("AGTTCGAATTATATGCT"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 7);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(g.getColorForSampleName("kid"))
                .displayColors(g.getColorsForSampleNames(Arrays.asList("mom", "dad", "kid")))
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.AND)
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
                .connectAllNeighbors(false)
                .stopper(ExplorationStopper.class)
                .graph(g)
                .make();

        String seed = "ATTATAT";

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> walk = e.dfs(seed);

        Set<CortexVertex> iss = new HashSet<>();
        Set<CortexVertex> oss = new HashSet<>();

        for (CortexVertex v : walk.vertexSet()) {
            if (TraversalEngine.inDegree(walk, v) > 1) {
                iss.add(v);
            }

            if (TraversalEngine.outDegree(walk, v) > 1) {
                oss.add(v);
            }
        }

        e.getConfiguration().setStoppingRule(DestinationStopper.class);
        e.getConfiguration().setGraphCombinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.OR);

        for (CortexVertex is : iss) {
            for (CortexVertex os : oss) {
                DirectedGraph<CortexVertex, CortexEdge> p = new DefaultDirectedGraph<>(CortexEdge.class);
                p.addVertex(os);

                e.getConfiguration().setPreviousTraversal(p);

                for (int c : g.getColorsForSampleNames(Arrays.asList("mom", "dad"))) {
                    e.getConfiguration().setTraversalColor(c);

                    DirectedWeightedPseudograph<CortexVertex, CortexEdge> w = e.dfs(is.getSk());

                    if (w != null) {
                        Graphs.addGraph(walk, w);
                    }
                }
            }
        }

        PathFinder dspKid = new PathFinder(walk, g.getColorForSampleName("kid"));
        PathFinder dspMom = new PathFinder(walk, g.getColorForSampleName("mom"));
        PathFinder dspDad = new PathFinder(walk, g.getColorForSampleName("dad"));

        Set<Pair<String, String>> allelesVsMom = new HashSet<>();
        Set<Pair<String, String>> allelesVsDad = new HashSet<>();

        for (CortexVertex os : oss) {
            for (CortexVertex is : iss) {
                GraphPath<CortexVertex, CortexEdge> pk = dspKid.getPathFinder(os, is, new CortexKmer(seed), true);
                GraphPath<CortexVertex, CortexEdge> pm = dspMom.getPathFinder(os, is);
                GraphPath<CortexVertex, CortexEdge> pd = dspDad.getPathFinder(os, is);

                Pair<String, String> akm = pathsToAlleles(pk, pm);
                Pair<String, String> akd = pathsToAlleles(pk, pd);

                allelesVsMom.add(akm);
                allelesVsDad.add(akd);
            }
        }

        Pair<String, String> alleleVsMom = allelesVsMom.iterator().next();
        Pair<String, String> alleleVsDad = allelesVsDad.iterator().next();

        Assert.assertEquals(alleleVsMom.getFirst(), "");
        Assert.assertEquals(alleleVsMom.getSecond(), "CTGGGC");
        Assert.assertEquals(alleleVsDad.getFirst(), "");
        Assert.assertEquals(alleleVsDad.getSecond(), "CTGGGC");
    }

    @Test
    public void testCallHeterozygote() {
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom",                Collections.singletonList("AGTTCGAATCTGGGCTATATGCT"));
        haplotypes.put("dad",                Collections.singletonList("AGTTCGAATCTGGGCTATATGCT"));
        haplotypes.put("kid", Arrays.asList("AGTTCGAATCTGGGCTATATGCT", "AGTTCGAATCTGAGCTATATGCT"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 7);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(g.getColorForSampleName("kid"))
                .displayColors(g.getColorsForSampleNames(Arrays.asList("mom", "dad", "kid")))
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.AND)
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
                .connectAllNeighbors(false)
                .stopper(ExplorationStopper.class)
                .graph(g)
                .make();

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> walk = e.dfs("TGAGCTA");

        Set<CortexVertex> iss = new HashSet<>();
        Set<CortexVertex> oss = new HashSet<>();

        for (CortexVertex v : walk.vertexSet()) {
            if (TraversalEngine.inDegree(walk, v) > 1) {
                iss.add(v);
            }

            if (TraversalEngine.outDegree(walk, v) > 1) {
                oss.add(v);
            }
        }

        e.getConfiguration().setStoppingRule(DestinationStopper.class);
        e.getConfiguration().setGraphCombinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.OR);

        for (CortexVertex is : iss) {
            for (CortexVertex os : oss) {
                DirectedGraph<CortexVertex, CortexEdge> p = new DefaultDirectedGraph<>(CortexEdge.class);
                p.addVertex(os);

                e.getConfiguration().setPreviousTraversal(p);

                for (int c : g.getColorsForSampleNames(Arrays.asList("mom", "dad"))) {
                    e.getConfiguration().setTraversalColor(c);

                    DirectedWeightedPseudograph<CortexVertex, CortexEdge> w = e.dfs(is.getSk());

                    if (w != null) {
                        Graphs.addGraph(walk, w);
                    }
                }
            }
        }

        PathFinder dspKid = new PathFinder(walk, g.getColorForSampleName("kid"));

        Set<Pair<String, String>> allelesVsSelf = new HashSet<>();

        for (CortexVertex os : oss) {
            for (CortexVertex is : iss) {
                GraphPath<CortexVertex, CortexEdge> pk1 = dspKid.getPathFinder(os, is, new CortexKmer("TGAGCTA"), true);
                GraphPath<CortexVertex, CortexEdge> pk2 = dspKid.getPathFinder(os, is, new CortexKmer("TGAGCTA"), false);

                Pair<String, String> akm = pathsToAlleles(pk1, pk2);
                allelesVsSelf.add(akm);
            }
        }

        Pair<String, String> alleleVsSelf = allelesVsSelf.iterator().next();

        Assert.assertEquals(alleleVsSelf.getFirst(), "A");
        Assert.assertEquals(alleleVsSelf.getSecond(), "G");
    }

    @Test
    public void testCallNAHR() {
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();

        String var1 = new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(30));
        String var2 = new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(30));
        String nahr = var1.substring(0, 15) + var2.substring(15, 30);

        System.out.println(var1);
        System.out.println(var2);
        System.out.println(nahr);

        haplotypes.put("var", Arrays.asList(var1, var2));
        haplotypes.put("kid", Collections.singletonList(nahr));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 13);

        TraversalEngine e = new TraversalEngineFactory()
                //.traversalColor(g.getColorForSampleName("kid"))
                .displayColors(g.getColorsForSampleNames(Arrays.asList("var", "kid")))
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.AND)
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
                .connectAllNeighbors(false)
                .stopper(ExplorationStopper.class)
                .graph(g)
                .make();

        String seed1 = var1.substring(7, 20);
        String seed2 = var2.substring(7, 20);
        String seedr = nahr.substring(7, 20);

        e.getConfiguration().setTraversalColor(g.getColorForSampleName("var"));
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> walk1 = e.dfs(seed1);
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> walk2 = e.dfs(seed2);

        e.getConfiguration().setTraversalColor(g.getColorForSampleName("kid"));
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> walkr = e.dfs(seedr);

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> graph = new DirectedWeightedPseudograph<>(CortexEdge.class);
        Graphs.addGraph(graph, walk1);
        Graphs.addGraph(graph, walk2);
        Graphs.addGraph(graph, walkr);

        GraphVisualizer gv = new GraphVisualizer(9000);
        gv.display(graph, "nahr");

        /*
        Set<CortexVertex> iss = new HashSet<>();
        Set<CortexVertex> oss = new HashSet<>();

        for (CortexVertex v : walk.vertexSet()) {
            if (TraversalEngine.inDegree(walk, v) > 1) {
                iss.add(v);
            }

            if (TraversalEngine.outDegree(walk, v) > 1) {
                oss.add(v);
            }
        }

        e.getConfiguration().setStoppingRule(DestinationStopper.class);
        e.getConfiguration().setGraphCombinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.OR);

        for (CortexVertex is : iss) {
            for (CortexVertex os : oss) {
                DirectedGraph<CortexVertex, CortexEdge> p = new DefaultDirectedGraph<>(CortexEdge.class);
                p.addVertex(os);

                e.getConfiguration().setPreviousTraversal(p);

                for (int c : g.getColorsForSampleNames(Arrays.asList("mom", "dad"))) {
                    e.getConfiguration().setTraversalColor(c);

                    DirectedWeightedPseudograph<CortexVertex, CortexEdge> w = e.dfs(is.getSk());

                    if (w != null) {
                        Graphs.addGraph(walk, w);
                    }
                }
            }
        }

        PathFinder dspKid = new PathFinder(walk, g.getColorForSampleName("kid"));

        Set<Pair<String, String>> allelesVsSelf = new HashSet<>();

        for (CortexVertex os : oss) {
            for (CortexVertex is : iss) {
                GraphPath<CortexVertex, CortexEdge> pk1 = dspKid.getPathFinder(os, is, new CortexKmer("TGAGCTA"), true);
                GraphPath<CortexVertex, CortexEdge> pk2 = dspKid.getPathFinder(os, is, new CortexKmer("TGAGCTA"), false);

                Pair<String, String> akm = pathsToAlleles(pk1, pk2);
                allelesVsSelf.add(akm);
            }
        }

        Pair<String, String> alleleVsSelf = allelesVsSelf.iterator().next();

        Assert.assertEquals(alleleVsSelf.getFirst(), "A");
        Assert.assertEquals(alleleVsSelf.getSecond(), "G");
        */
    }

    private static Pair<String, String> pathsToAlleles(GraphPath<CortexVertex, CortexEdge> p0, GraphPath<CortexVertex, CortexEdge> p1) {
        String s0 = pathToString(p0);
        String s1 = pathToString(p1);

        int s0start = 0, s0end = s0.length();
        int s1start = 0, s1end = s1.length();

        for (int i = 0, j = 0; i < s0.length() && j < s1.length(); i++, j++) {
            if (s0.charAt(i) != s1.charAt(j)) {
                s0start = i;
                s1start = j;
                break;
            }
        }

        for (int i = s0.length() - 1, j = s1.length() - 1; i >= 0 && j >= 0; i--, j--) {
            if (s0.charAt(i) != s1.charAt(j)) {
                s0end = i + 1;
                s1end = j + 1;
                break;
            }
        }

        return new Pair<>(s0.substring(s0start, s0end), s1.substring(s1start, s1end));
    }

    private static String pathToString(GraphPath<CortexVertex, CortexEdge> pk) {
        StringBuilder sbk = new StringBuilder();
        if (pk != null) {
            for (CortexVertex v : pk.getVertexList()) {
                if (sbk.length() == 0) {
                    sbk.append(v.getSk());
                } else {
                    sbk.append(v.getSk().substring(v.getSk().length() - 1, v.getSk().length()));
                }
            }
        }
        return sbk.toString();
    }
}
