package uk.ac.ox.well.indiana.commands.caller.call;

import org.jetbrains.annotations.NotNull;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.stoppingconditions.BubbleClosingStopper;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;

import java.io.File;
import java.io.PrintStream;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

/**
 * Created by kiran on 23/06/2017.
 */
public class Call extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="drafts", shortName="d", doc="Drafts")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Argument(fullName="reference", shortName="R", doc="Reference")
    public KmerLookup REFERENCE;

    @Argument(fullName="annotations", shortName="a", doc="Annotated contigs")
    public File ANNOTATIONS;

    @Output
    public PrintStream out;

    private Random rng = new Random();

    @Override
    public void execute() {
        //int childColor = GRAPH.getColorForSampleName(CHILD);
        //List<Integer> parentColors = GRAPH.getColorsForSampleNames(PARENTS);
        //List<Integer> recruitColors = GRAPH.getColorsForSampleNames(new ArrayList<>(LOOKUPS.keySet()));

        Map<String, List<Map<String, String>>> contigs = loadAnnotations();

        for (String contigName : contigs.keySet()) {
            log.info("Processing contig {}", contigName);

            callVariants(contigName, contigs.get(contigName));
        }
    }

    private void callVariants(String contigName, List<Map<String, String>> annotations) {
        for (String background : LOOKUPS.keySet()) {
            log.info("  background {}", background);

            String annotatedContig = annotateContig(annotations, background);

            int numTemplateSwitches = numTemplateSwitches(annotatedContig);
            int numNovelRuns = numNovelRuns(annotatedContig);

            numSimpleBubbles(annotatedContig, annotations);

            if (numTemplateSwitches >= 2 && numNovelRuns >= 2) {
                log.info("nahr: {} {} {} {}", contigName, numTemplateSwitches, numNovelRuns, annotatedContig);
            }
        }
    }

    private void numSimpleBubbles(String annotatedContig, List<Map<String, String>> annotations) {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(PARENTS);

        for (int c : parentColors) {
            log.info("    parent: {}", c);

            DirectedGraph<CortexVertex, CortexEdge> childContig = new DefaultDirectedGraph<>(CortexEdge.class);

            Set<String> traversalSeeds = new HashSet<>();

            CortexVertex cvl = null;
            //for (Map<String, String> ma : annotations) {

            String mostRecentNonNovelKmer = annotations.get(0).get("sk");
            boolean inNovelRun = false;

            for (int i = 0; i < annotations.size(); i++) {
                Map<String, String> ma = annotations.get(i);

                String sk = ma.get("sk");
                CortexKmer ck = new CortexKmer(sk);
                CortexRecord cr = GRAPH.findRecord(ck);

                CortexVertex cv = new CortexVertex(sk, cr);
                childContig.addVertex(cv);
                if (cvl != null) {
                    childContig.addEdge(cvl, cv, new CortexEdge(childColor, 1.0));
                }
                cvl = cv;

                log.info("before {} {} {}", i, inNovelRun, mostRecentNonNovelKmer);

                if (ma.get("is_novel").equals("false")) {
                    if (inNovelRun) {
                        traversalSeeds.add(sk);
                    }

                    mostRecentNonNovelKmer = sk;
                    inNovelRun = false;
                } else if (ma.get("is_novel").equals("true")) {
                    if (!inNovelRun) {
                        traversalSeeds.add(sk);
                    }

                    inNovelRun = true;
                }

                log.info("after {} {} {}", i, inNovelRun, mostRecentNonNovelKmer);

                /*
                Map<Integer, Set<String>> incomingKmers = TraversalEngine.getAllPrevKmers(cr, ck.isFlipped());
                Map<Integer, Set<String>> outgoingKmers = TraversalEngine.getAllNextKmers(cr, ck.isFlipped());

                incomingKmers.get(c).removeAll(incomingKmers.get(childColor));
                outgoingKmers.get(c).removeAll(outgoingKmers.get(childColor));

                traversalSeeds.addAll(incomingKmers.get(c));
                traversalSeeds.addAll(outgoingKmers.get(c));
                */
            }

            log.info("      num seeds: {}", traversalSeeds.size());

            TraversalEngine e = new TraversalEngineFactory()
                    .combinationOperator(OR)
                    .traversalDirection(BOTH)
                    .traversalColor(c)
                    .joiningColors(childColor)
                    .previousTraversal(childContig)
                    .stopper(BubbleClosingStopper.class)
                    .graph(GRAPH)
                    .make();

            DirectedWeightedPseudograph<CortexVertex, CortexEdge> gAll = new DirectedWeightedPseudograph<>(CortexEdge.class);
            for (String seed : traversalSeeds) {
                DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = e.dfs(seed);
                if (g != null && g.vertexSet().size() > 0) {
                    Graphs.addGraph(gAll, g);
                }
            }
        }
    }

    private int numTemplateSwitches(String annotatedContig) {
        final String flankingNovelRegex = "(([^_\\.])\\2+)_*(\\.+)_*(([^_\\.])\\5+)";
        Pattern flankingNovelPattern = Pattern.compile(flankingNovelRegex);

        Matcher flankingNovelMatcher = flankingNovelPattern.matcher(annotatedContig);
        int numSwitches = 0;
        if (flankingNovelMatcher.find()) {
            do {
                if (!flankingNovelMatcher.group(2).equals(flankingNovelMatcher.group(5))) {
                    numSwitches++;
                }
            } while (flankingNovelMatcher.find(flankingNovelMatcher.start(3)));
        }

        return numSwitches;
    }

    private int numNovelRuns(String annotatedContig) {
        final String novelRegex = "(\\.+)";
        Pattern novelPattern = Pattern.compile(novelRegex);

        Matcher novelMatcher = novelPattern.matcher(annotatedContig);
        int numNovelRuns = 0;
        while (novelMatcher.find()) {
            numNovelRuns++;
        }

        return numNovelRuns;
    }

    private String annotateContig(List<Map<String, String>> annotations, String background) {
        Map<String, String> chrCodes = createContigEncoding(annotations, background);

        StringBuilder ab = new StringBuilder();

        for (Map<String, String> m : annotations) {
            String[] lociStrings = m.get(background).split(";");

            String code = "_";

            if (m.get("is_novel").equals("true")) {
                code = ".";
            } else if (lociStrings.length == 1) {
                String[] pieces = lociStrings[0].split(":");
                String contig = pieces[0];

                code = chrCodes.get(contig);
            }

            ab.append(code);
        }

        return ab.toString();
    }

    @NotNull
    private Map<String, String> createContigEncoding(List<Map<String, String>> annotations, String background) {
        Map<String, String> contigEncoding = new HashMap<>();

        final String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
        Set<String> usedCodes = new HashSet<>();

        for (Map<String, String> m : annotations) {
            String[] lociStrings = m.get(background).split(";");

            if (lociStrings.length == 1) {
                String[] pieces = lociStrings[0].split(":");
                String contig = pieces[0];

                if (!contigEncoding.containsKey(contig)) {
                    String c;
                    do {
                        c = String.valueOf(alphabet.charAt(rng.nextInt(alphabet.length())));
                    } while(usedCodes.contains(c) && usedCodes.size() < alphabet.length());

                    contigEncoding.put(contig, c);
                    usedCodes.add(c);
                }
            }
        }

        return contigEncoding;
    }

    private Map<String, List<Map<String, String>>> loadAnnotations() {
        TableReader tr = new TableReader(ANNOTATIONS);

        Map<String, List<Map<String, String>>> contigs = new TreeMap<>();

        for (Map<String, String> m : tr) {
            if (!contigs.containsKey(m.get("name"))) {
                contigs.put(m.get("name"), new ArrayList<>());
            }
            contigs.get(m.get("name")).add(m);
        }

        return contigs;
    }
}
