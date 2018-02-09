package uk.ac.ox.well.cortexjdk.utils.traversal;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
import org.jgrapht.GraphPath;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.testng.Assert;
import org.testng.annotations.Test;
import uk.ac.ox.well.cortexjdk.CortexJDK;
import uk.ac.ox.well.cortexjdk.utils.assembler.TempGraphAssembler;
import uk.ac.ox.well.cortexjdk.utils.assembler.TempLinksAssembler;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinksIterable;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinksRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.CycleCollapsingContigStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.DestinationStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ExplorationStopper;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

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

        Map<String, String> expectations = new LinkedHashMap<>();
        expectations.put("mom", "AGTTCTGATCTGGGCTATATGCT");
        expectations.put("dad",   "TTCGAATCTGGGCTATATGCT");
        expectations.put("kid", "AGTTCTGATCTGGGCTATGGCT");

        for (int c = 0; c < 3; c++) {
            TraversalEngine e = new TraversalEngineFactory()
                    .traversalColor(c)
                    .graph(g)
                    .stoppingRule(ContigStopper.class)
                    .make();

            String contig = TraversalEngine.toContig(e.walk("CTGGG"));

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
                .traversalDirection(BOTH)
                .connectAllNeighbors(false)
                .stoppingRule(ContigStopper.class)
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

            String contig = TraversalEngine.toContig(e.walk("GTTCT"));

            Assert.assertEquals(contig, expectations.get(useRecruitment));
        }
    }

    @Test
    public void testCyclesWithoutLinksAreNotAssembled() {
        // This is the example in Figure 1 of the McCortex paper
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("test", Collections.singletonList("ACTGATTTCGATGCGATGCGATGCCACGGTGG"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 5);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(g.getColorForSampleName("test"))
                .stoppingRule(ContigStopper.class)
                .graph(g)
                .make();

        String contig = TraversalEngine.toContig(e.walk("ACTGA"));

        Assert.assertEquals(contig, "ACTGATTTCGATGC");
    }

    @Test
    public void testCyclesWithLinksAreAssembled() {
        // This is the example in Figure 1 of the McCortex paper
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("test", Collections.singletonList("ACTGATTTCGATGCGATGCGATGCCACGGTGG"));

        Map<String, Collection<String>> reads = new LinkedHashMap<>();
        reads.put("test", Collections.singletonList("TTTCGATGCGATGCGATGCCACG"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 5);
        CortexLinks l = TempLinksAssembler.buildLinks(g, reads, "test");

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(g.getColorForSampleName("test"))
                .stoppingRule(ContigStopper.class)
                .graph(g)
                .links(l)
                .make();

        String contig = TraversalEngine.toContig(e.walk("ACTGA"));

        Assert.assertEquals(contig, "ACTGATTTCGATGCGATGCGATGCCACGGTGG");
    }

    @Test
    public void iterateFwdWithoutPathInformation() {
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom",                Collections.singletonList("AGTTCGAATCTGGGCTATATGCT"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 7);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(g.getColorForSampleName("mom"))
                .graph(g)
                .make();

        String sk = "AGTTCGA";
        StringBuilder sb = new StringBuilder();
        sb.append(sk);

        e.seek(sk);
        while (e.hasNext()) {
            CortexVertex cv = e.next();
            sk = cv.getKmerAsString();

            sb.append(sk.substring(sk.length() - 1, sk.length()));
        }

        Assert.assertEquals(sb.toString(), haplotypes.get("mom").iterator().next());
    }

    @Test
    public void iterateRevWithoutPathInformation() {
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("mom",                Collections.singletonList("AGTTCGAATCTGGGCTATATGCT"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 7);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(g.getColorForSampleName("mom"))
                .graph(g)
                .make();

        String sk = "ATATGCT";
        StringBuilder sb = new StringBuilder();
        sb.append(sk);

        e.seek(sk);
        while (e.hasPrevious()) {
            CortexVertex cv = e.previous();
            sk = cv.getKmerAsString();

            sb.insert(0, sk.substring(0, 1));
        }

        Assert.assertEquals(sb.toString(), haplotypes.get("mom").iterator().next());
    }

    @Test
    public void iterateFwdToForkWithoutPathInformation() {
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("kid", Arrays.asList("AGTTCGAATCTGGGCTATATGCT", "AGTTCGAATCTGAGCTATATGCT"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 7);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(g.getColorForSampleName("kid"))
                .graph(g)
                .make();

        String sk = "AGTTCGA";
        StringBuilder sb = new StringBuilder();
        sb.append(sk);

        e.seek(sk);
        while (e.hasNext()) {
            CortexVertex cv = e.next();
            sk = cv.getKmerAsString();

            sb.append(sk.substring(sk.length() - 1, sk.length()));
        }

        Assert.assertEquals(sb.toString(), "AGTTCGAATCTG");
    }

    @Test
    public void iterateRevToForkWithoutPathInformation() {
        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("kid", Arrays.asList("AGTTCGAATCTGGGCTATATGCT", "AGTTCGAATCTGAGCTATATGCT"));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, 7);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(g.getColorForSampleName("kid"))
                .graph(g)
                .make();

        String sk = "ATATGCT";
        StringBuilder sb = new StringBuilder();
        sb.append(sk);

        e.seek(sk);
        while (e.hasPrevious()) {
            CortexVertex cv = e.previous();
            sk = cv.getKmerAsString();

            sb.insert(0, sk.substring(0, 1));
        }

        Assert.assertEquals(sb.toString(), "GCTATATGCT");
    }

    @Test
    public void testGoForwardAndBackward() {
        String hap = "AGTTCGAATCTGAGCTATATGCT";
        int kmerSize = 7;

        Map<String, Collection<String>> haplotypes = new LinkedHashMap<>();
        haplotypes.put("kid", Arrays.asList(hap));

        CortexGraph g = TempGraphAssembler.buildGraph(haplotypes, kmerSize);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(g.getColorForSampleName("kid"))
                .graph(g)
                .make();

        for (int i = 1; i < hap.length() - kmerSize; i++) {
            String sk = hap.substring(i, i + kmerSize);

            e.seek(sk);
            if (e.hasPrevious() && e.hasNext()) {
                e.next();
                CortexVertex cv = e.previous();

                Assert.assertEquals(cv.getKmerAsString(), sk);
            }
        }
    }

}
