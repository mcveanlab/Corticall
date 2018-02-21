package uk.ac.ox.well.cortexjdk.utils.traversal;

import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.junit.Test;
import org.testng.Assert;
import uk.ac.ox.well.cortexjdk.utils.assembler.TempGraphAssembler;
import uk.ac.ox.well.cortexjdk.utils.assembler.TempLinksAssembler;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;

import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

public class TraversalUtilsTest {
    @Test
    public void testGapFilling() {
        Map<String, Collection<String>> haplotypes = new HashMap<>();
        haplotypes.put("mom", Arrays.asList("TGGCTAGGTCATTATGATATTAAAATGCTAGCGC"));
        haplotypes.put("kid", Arrays.asList("TGGCTAGGTCATTATGAGATTAAAATGCTAGCGC"));

        CortexGraph g    = TempGraphAssembler.buildGraph(haplotypes, 7);
        CortexLinks lmom = TempLinksAssembler.buildLinks(g, haplotypes, "mom");
        CortexLinks lkid = TempLinksAssembler.buildLinks(g, haplotypes, "kid");

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(g.getColorForSampleName("kid"))
                .traversalDirection(BOTH)
                .combinationOperator(OR)
                .stoppingRule(ContigStopper.class)
                .graph(g)
                .links(lkid)
                .make();

        List<CortexVertex> w = e.walk("TGAGATT");

        Assert.assertEquals(TraversalUtils.toContig(w), haplotypes.get("kid").iterator().next());

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gGapped = TraversalUtils.toGraph(w);
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gFilled = TraversalUtils.fillGaps(gGapped, g, Arrays.asList(lmom, lkid));

        String contigKid = TraversalUtils.toContig(TraversalUtils.toWalk(gFilled, "TGGCTAG", g.getColorForSampleName("kid")));
        String contigMom = TraversalUtils.toContig(TraversalUtils.toWalk(gFilled, "TGGCTAG", g.getColorForSampleName("mom")));

        Assert.assertEquals(contigKid, haplotypes.get("kid").iterator().next());
        Assert.assertEquals(contigMom, haplotypes.get("mom").iterator().next());
    }
}
