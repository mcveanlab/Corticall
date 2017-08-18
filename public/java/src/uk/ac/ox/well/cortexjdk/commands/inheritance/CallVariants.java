package uk.ac.ox.well.cortexjdk.commands.inheritance;

import htsjdk.samtools.util.Interval;
import org.jetbrains.annotations.NotNull;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.*;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.stoppingconditions.BubbleOpeningStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingconditions.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexEdge;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

/**
 * Created by kiran on 08/08/2017.
 */
public class CallVariants extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Links")
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName="references", shortName="r", doc="References")
    public HashMap<String, KmerLookup> REFERENCES;

    @Argument(fullName="parent", shortName="p", doc="Parents")
    public HashMap<String, String> PARENTS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int refColor = GRAPH.getColorForSampleName("ref");
        Set<Integer> parentColors = new TreeSet<>(GRAPH.getColorsForSampleNames(new ArrayList<>(PARENTS.values())));
        Set<Integer> draftColors = new TreeSet<>(GRAPH.getColorsForSampleNames(new ArrayList<>(REFERENCES.keySet())));
        Set<Integer> childColors = getChildColors(parentColors, draftColors, refColor);

        log.info("Colors:");
        log.info("  - parents:  {}", parentColors);
        log.info("  - children: {}", childColors);
        log.info("  - refs:     {}", draftColors);

        Set<CortexKmer> seeds = getVariantSeeds(refColor, parentColors, draftColors);

        TraversalEngine e = new TraversalEngineFactory()
                .combinationOperator(AND)
                .traversalDirection(BOTH)
                .joiningColors(parentColors)
                .stopper(ContigStopper.class)
                .graph(GRAPH)
                .links(LINKS)
                .make();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Building contigs")
                .message("seeds processed")
                .maxRecord(seeds.size())
                .make(log);

        for (CortexKmer ck : seeds) {
            for (int c : childColors) {
                CortexRecord cr = GRAPH.findRecord(ck);
                if (cr.getCoverage(c) > 0) {
                    e.getConfiguration().setTraversalColor(c);

                    List<CortexVertex> cvs = e.walk(ck.getKmerAsString());

                    for (CortexVertex cv : cvs) {
                    }

                    break;
                }
            }

            pm.update();

            /*
            CortexRecord cr = seeds.get(ck);

            for (int c : childColors) {
                if (cr.getCoverage(c) > 0) {
                    Collection<String> inEdges = cr.getInEdgesAsStrings(c, false);
                    for (String inEdge : inEdges) {
                        CortexKmer inKmer = new CortexKmer(inEdge + cr.getKmerAsString().substring(0, GRAPH.getKmerSize() - 1));


                    }

                    Collection<String> outEdges = cr.getOutEdgesAsStrings(c, false);
                }
            }
            */
        }




                //sg.addVertex(new CortexVertex(new CortexByteKmer(cr.getKmerAsBytes()), cr));

                //canWalkToCanonicalReference(cr, childColors, refColor)) {

////                int draftColor = getDraftColor(cr, draftColors);
////                KmerLookup kl = REFERENCES.get(GRAPH.getSampleName(draftColor));
////
////                log.info("{} {}", cr.getKmerAsString(), kl.findKmer(cr.getKmerAsString()));
//
//                seeds++;
//
//                /*
//                List<Interval> intervals = getCanonicalReferenceCoordinates(cr, childColors, refColor);
//
//                if (intervals != null) {
//                    int draftColor = getDraftColor(cr, draftColors);
//                    int bubbleColor = getBubbleColor(draftColors, draftColor);
//
//                    Set<String> childAllele = new HashSet<>();
//                    Set<String> draftAllele = new HashSet<>();
//                    Set<Interval> locus = new HashSet<>();
//                    Set<Integer> colors = new HashSet<>();
//
//                    for (int cc : childColors) {
//                        if (cr.getCoverage(cc) > 0) {
//                            Pair<String, String> bubbleB = callSimpleBubble(cr, cc, bubbleColor);
//                            //Pair<String, String> bubbleD = callSimpleBubble(cr, cc, draftColor);
//
//                            if (bubbleB != null) {
//                                Interval coords = getBubbleCanonicalCoordinates(bubbleB.getFirst(), cc, refColor);
//
//                                if (coords != null) {
//                                    Pair<String, String> alleles = trimToAlleles(bubbleB);
//
//                                    //log.info("  - {} {} {} {}", cr.getKmerAsString(), alleles.getFirst(), alleles.getSecond(), coords);
//                                    if (alleles.getFirst().length() == alleles.getSecond().length() && alleles.getFirst().length() == 1) {
//                                        childAllele.add(alleles.getFirst());
//                                        draftAllele.add(alleles.getSecond());
//                                        locus.add(coords);
//                                        colors.add(cc);
//
//                                        for (int i = 0; i <= bubbleB.getFirst().length() - GRAPH.getKmerSize(); i++) {
//                                            String sk = bubbleB.getFirst().substring(i, i + GRAPH.getKmerSize());
//                                            seen.add(new CortexBinaryKmer(CortexRecord.encodeBinaryKmer(sk.getBytes())));
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                    }
//
//                    if (childAllele.size() == 1) {
//                        numVariants++;
//                        log.info("  call: {} {} {} {} {} {}", GRAPH.position(), GRAPH.getSampleName(draftColor), childAllele, draftAllele, locus, colors);
//
//                        Interval loc = locus.iterator().next();
//
//                        out.println(Joiner.on(" ").join(loc.getContig(), loc.getStart(), loc.getEnd(), draftColor, GRAPH.getSampleName(draftColor), childAllele, draftAllele, colors));
//                    }
//                }
//                */
//            }
//
//            pm.update("records processed, " + seeds + " seeds, " + numVariants + " variants, " + seen.size() + " kmers from variants");
//            pm.update();
//
//            pos++;
//            if (pos >= CHUNKSIZE) { break; }
//        }
    }

    @NotNull
    private Set<CortexKmer> getVariantSeeds(int refColor, Set<Integer> parentColors, Set<Integer> draftColors) {
        Set<CortexKmer> seeds = new HashSet<>();

        DirectedGraph<String, DefaultEdge> sg = new DefaultDirectedGraph<>(DefaultEdge.class);

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing records")
                .message("records processed")
                .maxRecord(GRAPH.getNumRecords())
                .updateRecord(GRAPH.getNumRecords() / 10)
                .make(log);

        for (CortexRecord cr : GRAPH) {
            if (isSinglyConnected(cr) &&
                isSharedWithOnlyOneParent(cr, parentColors, draftColors) &&
                isSharedWithSomeChildren(cr, parentColors, draftColors, refColor) &&
                hasUniqueCoordinates(cr, draftColors)) {

                seeds.add(cr.getCortexKmer());

                String skFwd = cr.getKmerAsString();
                String skRev = SequenceUtils.reverseComplement(cr.getKmerAsString());
                sg.addVertex(skFwd);
                sg.addVertex(skRev);

                for (int c = 0; c < cr.getNumColors(); c++) {
                    if (cr.getCoverage(c) > 0) {
                        Collection<String> ies = cr.getInEdgesAsStrings(c, false);
                        for (String ie : ies) {
                            String inEdgeFwd = ie + skFwd.substring(0, skFwd.length() - 1);
                            String outEdgeRev = SequenceUtils.reverseComplement(inEdgeFwd);

                            sg.addVertex(inEdgeFwd);
                            sg.addEdge(inEdgeFwd, skFwd);

                            sg.addVertex(outEdgeRev);
                            sg.addEdge(skRev, outEdgeRev);
                        }

                        Collection<String> oes = cr.getOutEdgesAsStrings(c, false);
                        for (String oe : oes) {
                            String outEdgeFwd = skFwd.substring(1, skFwd.length()) + oe;
                            String inEdgeRev = SequenceUtils.reverseComplement(outEdgeFwd);

                            sg.addVertex(outEdgeFwd);
                            sg.addEdge(skFwd, outEdgeFwd);

                            sg.addVertex(inEdgeRev);
                            sg.addEdge(inEdgeRev, skRev);
                        }
                    }
                }
            }

            pm.update();
        }

        Set<String> uniqueSeeds = new HashSet<>();
        Set<CortexKmer> goodSeeds = new HashSet<>();
        for (String sk : sg.vertexSet()) {
            if (sg.inDegreeOf(sk) == 0 && sg.outDegreeOf(sk) == 1) {
                uniqueSeeds.add(sk);

                List<String> contig = new ArrayList<>();
                contig.add(sk);

                String v = sk;
                while (sg.outDegreeOf(v) == 1) {
                    List<String> out = Graphs.successorListOf(sg, v);
                    v = out.get(0);

                    contig.add(v);
                }

                if (contig.size() > 3) {
                    goodSeeds.add(new CortexKmer(sk));
                }
            }
        }

        log.info("Found {} seed kmers for putative variants, {} unique, {} good", seeds.size(), uniqueSeeds.size(), goodSeeds.size());

        return goodSeeds;
    }

    /*
    private int getBubbleColor(Set<Integer> draftColors, int draftColor) {
        int bubbleColor = draftColor;
        for (int dc : draftColors) {
            if (dc != draftColor) {
                bubbleColor = dc;
            }
        }
        return bubbleColor;
    }
    */

    private Set<Integer> getChildColors(Set<Integer> parentColors, Set<Integer> draftColors, int refColor) {
        Set<Integer> childColors = new TreeSet<>();

        for (int c = 0; c < GRAPH.getNumColors(); c++) {
            childColors.add(c);
        }

        childColors.removeAll(parentColors);
        childColors.removeAll(draftColors);
        childColors.remove(refColor);

        return childColors;
    }

    private boolean isSharedWithOnlyOneParent(CortexRecord cr, Set<Integer> parentColors, Set<Integer> draftColors) {
        int numParentsWithCoverage = 0;
        for (int pc : parentColors) {
            if (cr.getCoverage(pc) > 0) { numParentsWithCoverage++; }
        }

        int numDraftsWithCoverage = 0;
        for (int dc : draftColors) {
            if (cr.getCoverage(dc) > 0) { numDraftsWithCoverage++; }
        }

        return numParentsWithCoverage == 1 && numDraftsWithCoverage == 1;
    }

    private boolean isSharedWithSomeChildren(CortexRecord cr, Set<Integer> parentColors, Set<Integer> draftColors, int refColor) {
        Set<Integer> ignoreColors = new HashSet<>();
        ignoreColors.addAll(parentColors);
        ignoreColors.addAll(draftColors);
        ignoreColors.add(refColor);

        int numChildrenWithCoverage = 0;
        int numChildren = 0;
        for (int c = 0; c < cr.getNumColors(); c++) {
            if (!ignoreColors.contains(c)) {
                if (cr.getCoverage(c) > 0) {
                    numChildrenWithCoverage++;
                }
                numChildren++;
            }
        }

        return numChildrenWithCoverage > 1; // && numChildrenWithCoverage < numChildren;
    }

    private boolean isSinglyConnected(CortexRecord cr) {
        boolean isSinglyConnected = true;
        for (int c = 0; c < cr.getNumColors(); c++) {
            if (cr.getCoverage(c) > 0) {
                if (!(cr.getInDegree(c) == 1 && cr.getOutDegree(c) == 1)) {
                    isSinglyConnected = false;
                }
            }
        }

        return isSinglyConnected;
    }

    private boolean hasUniqueCoordinates(CortexRecord cr, Set<Integer> draftColors) {
        int draftColor = -1;
        int numDraftsWithCoverage = 0;
        for (int dc : draftColors) {
            if (cr.getCoverage(dc) > 0) {
                draftColor = dc;
                numDraftsWithCoverage++;
            }
        }

        if (numDraftsWithCoverage == 1 && draftColor > -1) {
            KmerLookup kl = REFERENCES.get(GRAPH.getSampleName(draftColor));

            Set<Interval> its = kl.findKmer(cr.getKmerAsString());

            return its != null && its.size() == 1;
        }

        return false;
    }

    private boolean canWalkToCanonicalReference(CortexRecord cr, Set<Integer> childColors, int refColor) {
        TraversalEngine e = new TraversalEngineFactory()
                .graph(GRAPH)
                .make();

        for (int cc : childColors) {
            e.getConfiguration().setTraversalColor(cc);

            int found = 0;

            for (boolean goForward : Arrays.asList(true, false)) {
                e.seek(cr.getKmerAsString());
                while (goForward ? e.hasNext() : e.hasPrevious()) {
                    CortexVertex cv = goForward ? e.next() : e.previous();

                    if (cv.getCr().getCoverage(refColor) > 0) {
                        found++;
                        break;
                    }
                }
            }

            if (found == 2) {
                return true;
            }
        }

        return false;
    }

/*
    private List<Interval> getCanonicalReferenceCoordinates(CortexRecord cr, Set<Integer> childColors, int refColor) {
        TraversalEngine e = new TraversalEngineFactory()
                .graph(GRAPH)
                .make();

        for (int cc : childColors) {
            e.getConfiguration().setTraversalColor(cc);

            List<Interval> intervals = new ArrayList<>();

            for (boolean goForward : Arrays.asList(true, false)) {
                e.seek(cr.getKmerAsString());
                while (goForward ? e.hasNext() : e.hasPrevious()) {
                    CortexVertex cv = goForward ? e.next() : e.previous();

                    if (cv.getCr().getCoverage(refColor) > 0) {
                        Set<Interval> its = REFERENCE.findKmer(cv.getSk());
                        if (its != null && its.size() == 1) {
                            intervals.add(its.iterator().next());
                        }

                        break;
                    }
                }
            }

            if (intervals.size() == 2) {
                Interval it0 = intervals.get(0);
                Interval it1 = intervals.get(1);

                if (it0.getContig().equals(it1.getContig())) {
                    return intervals;
                }
            }
        }

        return null;
    }

    private int getDraftColor(CortexRecord cr, Set<Integer> draftColors) {
        for (int dc : draftColors) {
            if (cr.getCoverage(dc) > 0) {
                return dc;
            }
        }

        return -1;
    }

    private Pair<String, String> callSimpleBubble(CortexRecord cr, int childColor, int bubbleColor) {
        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(childColor)
                .graph(GRAPH)
                .make();

        List<CortexVertex> childVertices = new ArrayList<>();
        childVertices.add(new CortexVertex(cr.getKmerAsByteKmer(), cr));

        for (boolean goForward : Arrays.asList(false, true)) {
            e.seek(cr.getKmerAsString());
            while (goForward ? e.hasNext() : e.hasPrevious()) {
                CortexVertex cv = goForward ? e.next() : e.previous();

                if (goForward) { childVertices.add(cv); }
                else { childVertices.add(0, cv); }

                if (cv.getCr().getCoverage(bubbleColor) > 0) { break; }
            }
        }

        if (childVertices.get(0).getCr().getCoverage(bubbleColor) > 0) {
            e = new TraversalEngineFactory()
                    .traversalColor(bubbleColor)
                    .graph(GRAPH)
                    .make();

            List<CortexVertex> draftVertices = new ArrayList<>();
            draftVertices.add(childVertices.get(0));

            e.seek(childVertices.get(0).getSk());
            while (e.hasNext()) {
                CortexVertex cv = e.next();

                draftVertices.add(cv);

                if (cv.getCr().getCoverage(childColor) > 0) { break; }
            }

            if (draftVertices.get(draftVertices.size() - 1).getSk().equals(childVertices.get(childVertices.size() - 1).getSk())) {
                return new Pair<>(TraversalEngine.toContig(childVertices), TraversalEngine.toContig(draftVertices));
            }
        }

        return null;
    }

    private Pair<String, String> trimToAlleles(Pair<String, String> haplotypes) {
        String s0 = haplotypes.getFirst();
        String s1 = haplotypes.getSecond();

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
            if (s0.charAt(i) != s1.charAt(j) || i == s0start - 1 || j == s1start - 1) {
                s0end = i + 1;
                s1end = j + 1;
                break;
            }
        }

        String[] pieces = new String[4];
        pieces[0] = s0.substring(0, s0start);
        pieces[1] = s0.substring(s0start, s0end);
        pieces[2] = s1.substring(s1start, s1end);
        pieces[3] = s0.substring(s0end, s0.length());

        return new Pair<>(pieces[1], pieces[2]);
    }

    private Interval getBubbleCanonicalCoordinates(String childHaplotype, int childColor, int refColor) {
        Interval is = null, ie = null;

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(childColor)
                .graph(GRAPH)
                .make();

        String firstKmer = childHaplotype.substring(0, GRAPH.getKmerSize());
        String lastKmer = childHaplotype.substring(childHaplotype.length() - GRAPH.getKmerSize(), childHaplotype.length());

        for (boolean goForward : Arrays.asList(false, true)) {
            e.seek(goForward ? lastKmer : firstKmer);
            while (goForward ? e.hasNext() : e.hasPrevious()) {
                CortexVertex cv = goForward ? e.next() : e.previous();

                if (cv.getCr().getCoverage(refColor) > 0) {
                    Set<Interval> intervals = REFERENCE.findKmer(cv.getSk());

                    if (intervals != null && intervals.size() == 1) {
                        if (goForward) {
                            ie = intervals.iterator().next();
                        } else {
                            is = intervals.iterator().next();
                        }

                        break;
                    }
                }
            }
        }

        if (is != null && ie != null) {
            return is.getStart() < ie.getStart() ? is : ie;
        }

        return null;
    }
    */
}
