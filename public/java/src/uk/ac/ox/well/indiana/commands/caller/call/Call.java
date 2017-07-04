package uk.ac.ox.well.indiana.commands.caller.call;

import com.google.common.base.Joiner;
import org.jetbrains.annotations.NotNull;
import org.jgrapht.DirectedGraph;
import org.jgrapht.GraphPath;
import org.jgrapht.alg.shortestpath.DijkstraShortestPath;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
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
            if (isNahrEvent(contigName, contigs.get(contigName))) {
                //out.println(contigName);
            }
        }
    }

    private boolean isNahrEvent(String contigName, List<Map<String, String>> annotations) {
        String annotatedContig = annotateContig(annotations);

        int numTemplateSwitches = numTemplateSwitches(annotatedContig);
        int numNovelRuns = numNovelRuns(annotatedContig);

        //numBubbles(annotatedContig, annotations);

        log.info("{} {} {} {}", contigName, numTemplateSwitches, numNovelRuns, annotatedContig.length());

        if (numTemplateSwitches >= 2 && numNovelRuns >= 2) {
            log.info("nahr: {} {} {} {}", contigName, numTemplateSwitches, numNovelRuns, annotatedContig);

            out.println(contigName + "\t" + annotatedContig);

            return true;
        }

        return false;
    }

    private void numBubbles(String annotatedContig, List<Map<String, String>> annotations) {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(PARENTS);

        for (int c : parentColors) {
            log.info("    parent: {}", c);

            DirectedGraph<CortexVertex, CortexEdge> childContig = new DefaultDirectedGraph<>(CortexEdge.class);

            Set<CortexVertex> traversalSeeds = new HashSet<>();

            CortexVertex cvl = null;
            //for (Map<String, String> ma : annotations) {

            CortexVertex mostRecentNonNovelKmer = null;
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

                if (mostRecentNonNovelKmer == null) {
                    mostRecentNonNovelKmer = cv;
                }

                if (ma.get("is_novel").equals("false")) {
                    if (inNovelRun) {
                        traversalSeeds.add(cv);
                    }

                    mostRecentNonNovelKmer = cv;
                    inNovelRun = false;
                } else if (ma.get("is_novel").equals("true")) {
                    if (!inNovelRun) {
                        traversalSeeds.add(mostRecentNonNovelKmer);
                    }

                    inNovelRun = true;
                }
            }

            log.info("      num seeds: {}", traversalSeeds.size());

            TraversalEngine e = new TraversalEngineFactory()
                    .combinationOperator(OR)
                    .traversalDirection(BOTH)
                    .traversalColor(c)
                    .joiningColors(childColor)
                    //.previousTraversal(childContig)
                    .stopper(BubbleClosingStopper.class)
                    .graph(GRAPH)
                    .make();

            DirectedGraph<CortexVertex, CortexEdge> sg = new DefaultDirectedGraph<>(CortexEdge.class);
            for (CortexVertex seed : traversalSeeds) {
                sg.addVertex(seed);
            }

            for (CortexVertex source : traversalSeeds) {
                for (CortexVertex sink : traversalSeeds) {
                    if (!source.equals(sink)) {
                        e.getConfiguration().setPreviousTraversal(sg);

                        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = e.dfs(source.getSk());
                        if (g != null && g.vertexSet().size() > 0) {
                            DijkstraShortestPath<CortexVertex, CortexEdge> dsp = new DijkstraShortestPath<>(g);

                            if (!source.equals(sink) && g.containsVertex(source) && g.containsVertex(sink)) {
                                GraphPath<CortexVertex, CortexEdge> p = dsp.getPath(source, sink);

                                if (p != null) {
                                    log.info("{} {} {} {} {}", source, sink, p.getLength(), p.getStartVertex(), p.getEndVertex());
                                }
                            }
                        }
                    }
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

    private String annotateContig(List<Map<String, String>> annotations) {
        List<List<String>> annotatedContigs = new ArrayList<>();
        Map<String, String> chrCodes = createContigEncoding(annotations, LOOKUPS.keySet());

        for (String background : LOOKUPS.keySet()) {
            String annotatedContig = annotateContig(annotations, background, chrCodes);

            log.debug("{}", annotatedContig);

            String[] pieces = annotatedContig.split("((?<=\\.+)|(?=\\.+))");

            annotatedContigs.add(Arrays.asList(pieces));
        }

        List<String> finalPieces = new ArrayList<>();
        int lastBestBackgroundIndex = 0;
        for (int fragmentIndex = 0; fragmentIndex < annotatedContigs.get(0).size(); fragmentIndex++) {
            if (annotatedContigs.get(0).get(fragmentIndex).contains(".")) {
                finalPieces.add(annotatedContigs.get(0).get(fragmentIndex));
            } else {
                int bestBackgroundIndex = 0;
                int numKmersUniquelyPlaced = 0;

                for (int backgroundIndex = 0; backgroundIndex < annotatedContigs.size(); backgroundIndex++) {
                    Map<String, Integer> codeUsageMap = new HashMap<>();

                    for (int codeIndex = 0; codeIndex < annotatedContigs.get(backgroundIndex).get(fragmentIndex).length(); codeIndex++) {
                        char code = annotatedContigs.get(backgroundIndex).get(fragmentIndex).charAt(codeIndex);

                        if (code != '.' && code != '_') {
                            ContainerUtils.increment(codeUsageMap, String.valueOf(code));
                        }
                    }

                    String mostCommonCode = ContainerUtils.mostCommonKey(codeUsageMap);
                    if (mostCommonCode != null) {
                        int codeCount = codeUsageMap.get(mostCommonCode);

                        if (codeCount == numKmersUniquelyPlaced) {
                            bestBackgroundIndex = lastBestBackgroundIndex;
                        } else if (codeCount > numKmersUniquelyPlaced) {
                            bestBackgroundIndex = backgroundIndex;
                            numKmersUniquelyPlaced = codeCount;
                        }
                    } else {
                        // TODO: this is basically saying that within a single fragment, we weren't able to find any useful location information.  But other fragments may hold it.  Stopping here without looking ahead is dumb.  Fix this.
                        bestBackgroundIndex = 0;
                    }
                }

                finalPieces.add(annotatedContigs.get(bestBackgroundIndex).get(fragmentIndex));
                lastBestBackgroundIndex = bestBackgroundIndex;
            }
        }

        return Joiner.on("").join(finalPieces);
    }

    private String annotateContig(List<Map<String, String>> annotations, String background, Map<String, String> chrCodes) {
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
    private Map<String, String> createContigEncoding(List<Map<String, String>> annotations, Set<String> backgrounds) {
        Map<String, String> contigEncoding = new HashMap<>();

        final String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
        Set<String> usedCodes = new HashSet<>();

        for (Map<String, String> m : annotations) {
            for (String background : backgrounds) {
                String[] lociStrings = m.get(background).split(";");

                if (lociStrings.length == 1) {
                    String[] pieces = lociStrings[0].split(":");
                    String contig = pieces[0];

                    if (!contigEncoding.containsKey(contig)) {
                        String c;
                        do {
                            c = String.valueOf(alphabet.charAt(rng.nextInt(alphabet.length())));
                        } while (usedCodes.contains(c) && usedCodes.size() < alphabet.length());

                        contigEncoding.put(contig, c);
                        usedCodes.add(c);
                    }
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
