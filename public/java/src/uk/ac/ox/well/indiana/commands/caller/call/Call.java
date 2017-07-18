package uk.ac.ox.well.indiana.commands.caller.call;

import com.google.common.base.Joiner;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
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
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
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

        Map<String, List<Map<String, String>>> allAnnotations = loadAnnotations();
        for (String contigName : allAnnotations.keySet()) {
            String contig = getContig(allAnnotations.get(contigName));
            List<KmerAnnotation> annotatedContig = annotateContig(allAnnotations.get(contigName));
            List<KmerAnnotation> smoothedAnnotatedContig = smoothAnnotations(annotatedContig);

            log.info("{} {}", contigName, annotatedContig.size());

            for (int i = 0; i < annotatedContig.size(); i++) {
                log.info("  - {} {}", annotatedContig.get(i), smoothedAnnotatedContig.get(i));
            }

            //if (isNahrEvent(annotatedContig)) {
                //out.println(contigName);
            //}
        }
    }

    private boolean isNahrEvent(String annotatedContig) {
        int numTemplateSwitches = numTemplateSwitches(annotatedContig);
        int numNovelRuns = numNovelRuns(annotatedContig);

        if (numTemplateSwitches >= 1 && numNovelRuns >= 1) {
            return true;
        }

        return false;
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

    private String getContig(List<Map<String, String>> annotations) {
        StringBuilder sb = new StringBuilder();

        for (Map<String, String> m : annotations) {
            if (sb.length() == 0) {
                sb.append(m.get("sk"));
            } else {
                sb.append(m.get("sk").substring(GRAPH.getKmerSize() - 1, GRAPH.getKmerSize()));
            }
        }

        return sb.toString();
    }

    private class KmerAnnotation {
        private char code;
        private String background;
        private String intervals;
        private int offset;

        public KmerAnnotation() {}

        public KmerAnnotation(KmerAnnotation o) {
            code = o.getCode();
            background = o.getBackground();
            intervals = o.getIntervals();
            offset = o.getOffset();
        }

        public char getCode() { return code; }
        public void setCode(char code) { this.code = code; }

        public String getBackground() { return background; }
        public void setBackground(String background) { this.background = background; }

        public String getIntervals() { return intervals; }
        public void setIntervals(String intervals) { this.intervals = intervals; }

        public int getOffset() { return offset; }
        public void setOffset(int offset) { this.offset = offset; }

        @Override
        public String toString() {
            return "KmerAnnotation{" +
                    "code=" + code +
                    ", background='" + background + '\'' +
                    ", intervals='" + intervals + '\'' +
                    '}';
        }
    }

    private List<KmerAnnotation> smoothAnnotations(List<KmerAnnotation> annotatedContig) {
        List<KmerAnnotation> smoothedAnnotatedContig = new ArrayList<>();

        int lastValidIndex = -1;
        char lastValidCode = ' ';
        int nextValidIndex = -2;
        char nextValidCode = ' ';

        for (int i = 0; i < annotatedContig.size(); i++) {
            if (annotatedContig.get(i).getCode() == '_') {
                for (int j = i - 1; j >= 0; j--) {
                    char currentCode = annotatedContig.get(j).getCode();
                    if (currentCode != '_' && currentCode != '.' && currentCode != '?') {
                        lastValidIndex = j;
                        lastValidCode = currentCode;
                        break;
                    }
                }

                for (int j = i + 1; j < annotatedContig.size(); j++) {
                    char currentCode = annotatedContig.get(j).getCode();
                    if (currentCode != '_' && currentCode != '.' && currentCode != '?') {
                        nextValidIndex = j;
                        nextValidCode = currentCode;
                        break;
                    }
                }

                if (lastValidCode == ' ' && nextValidCode != ' ') {

                } else if (lastValidCode != ' ' && nextValidCode == ' ') {

                } else if (lastValidCode == nextValidCode) {
                    for (int j = i; j < nextValidIndex; j++) {
                        KmerAnnotation ka = new KmerAnnotation(annotatedContig.get(j));
                        ka.setCode(lastValidCode);
                    }
                    i = nextValidIndex;
                }
            } else {
                KmerAnnotation ka = new KmerAnnotation(annotatedContig.get(i));
                smoothedAnnotatedContig.add(ka);
            }
        }


        return smoothedAnnotatedContig;
    }

    private List<KmerAnnotation> annotateContig(List<Map<String, String>> annotations) {
        List<List<String>> annotatedContigs = new ArrayList<>();
        Set<String> usedAlphabet = new HashSet<>();
        Map<Character, String> annotationsToBackground = new HashMap<>();

        for (String background : LOOKUPS.keySet()) {
            String annotatedContig = annotateContig(annotations, background, usedAlphabet);

            for (int i = 0; i < annotatedContig.length(); i++) {
                char ann = annotatedContig.charAt(i);
                if (ann != '.' && ann != '_' && ann != '?' && !annotationsToBackground.containsKey(ann)) {
                    annotationsToBackground.put(ann, background);
                }
            }

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

                        if (code != '.' && code != '_' && code != '?') {
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

        List<KmerAnnotation> kmerAnnotations = new ArrayList<>();

        int offset = 0;
        for (String finalPiece : finalPieces) {
            for (int i = 0; i < finalPiece.length(); i++, offset++) {
                char code = finalPiece.charAt(i);
                String background = annotationsToBackground.containsKey(code) ? annotationsToBackground.get(code) : null;

                KmerAnnotation ka = new KmerAnnotation();
                ka.setCode(code);
                ka.setBackground(background);
                ka.setIntervals(annotations.get(offset).get(background));
                ka.setOffset(offset);
                kmerAnnotations.add(ka);
            }
        }

        return kmerAnnotations;
    }

    private String annotateContig(List<Map<String, String>> annotations, String background, Set<String> usedAlphabet) {
        StringBuilder ab = new StringBuilder();

        final String alphabet = "!\"#$%&'()*+,-/0123456789:;<=>@ABCDEFGHIJKLMNOPQRSTUVWXYZ[]^`abcdefghijklmnopqrstuvwxyz{|}~";

        IntervalTreeMap<String> itm = new IntervalTreeMap<>();
        for (Map<String, String> m : annotations) {
            String[] lociStrings = m.get(background).split(";");

            String code = "_";

            if (m.get("is_novel").equals("true")) {
                code = ".";
            } else if (lociStrings.length == 1) {
                String[] pieces = lociStrings[0].split("[:-]");
                String contig = pieces[0];

                if (contig.equals("NA")) {
                    //code = String.valueOf(alphabet.charAt(rng.nextInt(alphabet.length())));
                    code = "?";
                } else {
                    int start = Integer.valueOf(pieces[1]);
                    int end = Integer.valueOf(pieces[2]);

                    Interval it = new Interval(contig, start, end);
                    Interval lit = new Interval(contig, start - 500, end + 500);

                    if (itm.containsOverlapping(lit)) {
                        code = itm.getOverlapping(lit).iterator().next();
                    } else {
                        do {
                            code = String.valueOf(alphabet.charAt(rng.nextInt(alphabet.length())));
                        } while (usedAlphabet.contains(code) && usedAlphabet.size() < alphabet.length());
                    }

                    if (usedAlphabet.size() >= alphabet.length()) {
                        throw new IndianaException("Out of codes to assign.");
                    }

                    itm.put(it, code);
                }
            }

            ab.append(code);
        }

        return ab.toString();
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
