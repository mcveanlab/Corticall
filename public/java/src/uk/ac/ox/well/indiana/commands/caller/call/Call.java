package uk.ac.ox.well.indiana.commands.caller.call;

import com.google.common.base.Joiner;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import org.apache.commons.math3.util.Pair;
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
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(PARENTS);
        //List<Integer> recruitColors = GRAPH.getColorsForSampleNames(new ArrayList<>(LOOKUPS.keySet()));

        out.println(Joiner.on("\t").join("contigName", "index", "kmer", "code", "background", "intervals", "is_smoothed", "altInDegree", "altOutDegree"));

        Map<String, List<Map<String, String>>> allAnnotations = loadAnnotations();
        for (String contigName : allAnnotations.keySet()) {
            String contig = getContig(allAnnotations.get(contigName));

            List<KmerAnnotation> annotatedContig = annotateContig(allAnnotations.get(contigName));
            List<KmerAnnotation> smoothedAnnotatedContig = smoothAnnotations(annotatedContig, allAnnotations.get(contigName));
            addEdgeAnnotations(childColor, parentColors, contig, smoothedAnnotatedContig);

            for (int i = 0; i < smoothedAnnotatedContig.size(); i++) {
                KmerAnnotation ka = smoothedAnnotatedContig.get(i);
                KmerAnnotation kb = annotatedContig.get(i);
                out.println(Joiner.on("\t").join(contigName, i, ka.getKmer(), ka.getCode(), ka.getBackground(), ka.getIntervals(), ka.isSmoothed(), ka.getAltInDegree(), ka.getAltOutDegree(),
                                                 contigName, i, kb.getKmer(), kb.getCode(), kb.getBackground(), kb.getIntervals(), kb.isSmoothed(), kb.getAltInDegree(), kb.getAltOutDegree()));
            }

            //log.info("{} {} {} {}", contigName, contig.length(), contig.length() - GRAPH.getKmerSize() + 1, smoothedAnnotatedContig.size());

            //if (contigName.equals("contig291")) { }
        }
    }

    private void addEdgeAnnotations(int childColor, List<Integer> parentColors, String contig, List<KmerAnnotation> smoothedAnnotatedContig) {
        for (int i = 0; i < smoothedAnnotatedContig.size(); i++) {
            String sk = contig.substring(i, i + GRAPH.getKmerSize());
            CortexKmer ck = new CortexKmer(sk);
            CortexRecord cr = GRAPH.findRecord(ck);

            Map<Integer, Set<String>> inEdges = TraversalEngine.getAllPrevKmers(cr, ck.isFlipped());
            Map<Integer, Set<String>> outEdges = TraversalEngine.getAllNextKmers(cr, ck.isFlipped());

            smoothedAnnotatedContig.get(i).setKmer(sk);
            smoothedAnnotatedContig.get(i).setIncomingEdges(inEdges);
            smoothedAnnotatedContig.get(i).setOutgoingEdges(outEdges);

            for (int parentColor : parentColors) {
                Set<String> childIncomingEdges = new HashSet<>();
                childIncomingEdges.addAll(inEdges.get(childColor));
                childIncomingEdges.removeAll(inEdges.get(parentColor));

                if (childIncomingEdges.size() > 0) {
                    smoothedAnnotatedContig.get(i).setAltInDegree(childIncomingEdges.size());
                }

                Set<String> childOutgoingEdges = new HashSet<>();
                childOutgoingEdges.addAll(outEdges.get(childColor));
                childOutgoingEdges.removeAll(outEdges.get(parentColor));

                if (childOutgoingEdges.size() > 0) {
                    smoothedAnnotatedContig.get(i).setAltOutDegree(childOutgoingEdges.size());
                }
            }

            log.info("  {}", smoothedAnnotatedContig.get(i));
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
        private String kmer;
        private char code;
        private String background;
        private String intervals;
        private int offset;
        private boolean smoothed = false;
        private Map<Integer, Set<String>> incomingEdges = new HashMap<>();
        private Map<Integer, Set<String>> outgoingEdges = new HashMap<>();
        private int altInDegree = 0;
        private int altOutDegree = 0;

        public KmerAnnotation() {}

        public KmerAnnotation(KmerAnnotation o) {
            kmer = o.getKmer();
            code = o.getCode();
            background = o.getBackground();
            intervals = o.getIntervals();
            offset = o.getOffset();
            smoothed = o.isSmoothed();
            incomingEdges = o.getIncomingEdges();
            outgoingEdges = o.getOutgoingEdges();
        }

        public String getKmer() { return kmer == null ? "NA" : kmer; }
        public void setKmer(String kmer) { this.kmer = kmer; }

        public char getCode() { return code; }
        public void setCode(char code) { this.code = code; }

        public String getBackground() { return background == null ? "NA" : background; }
        public void setBackground(String background) { this.background = background; }

        public String getIntervals() { return intervals == null ? "NA" : intervals; }
        public void setIntervals(String intervals) { this.intervals = intervals; }

        public int getOffset() { return offset; }
        public void setOffset(int offset) { this.offset = offset; }

        public boolean isSmoothed() { return smoothed; }
        public void setSmoothed(boolean smoothed) { this.smoothed = smoothed; }

        public Map<Integer, Set<String>> getIncomingEdges() { return incomingEdges; }
        public void setIncomingEdges(Map<Integer, Set<String>> incomingEdges) { this.incomingEdges = incomingEdges; }

        public Map<Integer, Set<String>> getOutgoingEdges() { return outgoingEdges; }
        public void setOutgoingEdges(Map<Integer, Set<String>> outgoingEdges) { this.outgoingEdges = outgoingEdges; }

        public int getAltInDegree() { return altInDegree; }
        public void setAltInDegree(int altInDegree) { this.altInDegree = altInDegree; }

        public int getAltOutDegree() { return altOutDegree; }
        public void setAltOutDegree(int altOutDegree) { this.altOutDegree = altOutDegree; }

        @Override
        public String toString() {
            return "KmerAnnotation{" +
                    "kmer='" + kmer + '\'' +
                    ", code=" + code +
                    ", background='" + background + '\'' +
                    ", intervals='" + intervals + '\'' +
                    ", smoothed=" + smoothed +
                    ", i/o=" + altInDegree + "," + altOutDegree +
                    '}';
        }
    }

    private Interval stringToInterval(String intervalStr) {
        String[] pieces = intervalStr.split(":");
        String[] pos = pieces[1].split("-");

        return new Interval(pieces[0], Integer.valueOf(pos[0]), Integer.valueOf(pos[1]), pieces[2].equals("-"), "");
    }

    private String selectInterval(String manyIntervals, String singleInterval) {
        if (manyIntervals != null && singleInterval != null && !manyIntervals.contains("NA") && !singleInterval.contains("NA")) {
            Interval it = stringToInterval(singleInterval.split(";")[0]);
            Interval wideIt = new Interval(it.getContig(), it.getStart() - 500, it.getEnd() + 500, it.isNegativeStrand(), it.getName());

            for (String newIntervalStr : manyIntervals.split(";")) {
                Interval newInterval = stringToInterval(newIntervalStr);

                if (wideIt.intersects(newInterval)) {
                    return newIntervalStr;
                }
            }
        }

        return null;
    }

    private List<KmerAnnotation> smoothAnnotations(List<KmerAnnotation> annotatedContig, List<Map<String, String>> annotations) {
        List<KmerAnnotation> smoothedAnnotatedContig = new ArrayList<>();

        List<Pair<Integer, Integer>> lps = new ArrayList<>();
        int regionStart = -1;
        int regionEnd = -1;
        for (int i = 0; i < annotatedContig.size(); i++) {
            smoothedAnnotatedContig.add(new KmerAnnotation(annotatedContig.get(i)));

            if (annotatedContig.get(i).getCode() == '_') {
                if (regionStart == -1) {
                    regionStart = i;
                    regionEnd = i;
                } else if (regionEnd == i - 1) {
                    regionEnd = i;
                } else {
                    lps.add(new Pair<>(regionStart, regionEnd));
                    regionStart = i;
                    regionEnd = i;
                }
            }
        }

        if (regionStart > -1) {
            lps.add(new Pair<>(regionStart, regionEnd));
        }

        for (int i = 0; i < lps.size(); i++) {
            Pair<Integer, Integer> lp = lps.get(i);
            int lowerPosLimit = i > 0 ? lps.get(i-1).getSecond() + 1 : 0;
            int upperPosLimit = i < lps.size() - 1 ? lps.get(i+1).getFirst() - 1 : annotatedContig.size() - 1;
            int lastValidIndex = -1, nextValidIndex = -1;
            char lastValidCode = '_', nextValidCode = '_';

            for (int j = lp.getFirst() - 1; j >= lowerPosLimit; j--) {
                char currentCode = annotatedContig.get(j).getCode();
                if (currentCode == '.' || currentCode == '?') {
                    break;
                } else if (currentCode != '_') {
                    lastValidIndex = j;
                    lastValidCode = currentCode;
                    break;
                }
            }

            for (int j = lp.getSecond() + 1; j <= upperPosLimit; j++) {
                char currentCode = annotatedContig.get(j).getCode();
                if (currentCode == '.' || currentCode == '?') {
                    break;
                } else if (currentCode != '_') {
                    nextValidIndex = j;
                    nextValidCode = currentCode;
                    break;
                }
            }

            char fillCode = '_';
            int fillIndex = -1;
            if (lastValidIndex >= 0 && nextValidIndex >= 0 && lastValidCode == nextValidCode && lastValidCode != '_') {
                fillCode = lastValidCode;
                fillIndex = lastValidIndex;
            } else if (lastValidIndex == -1 && nextValidIndex >= 0 && nextValidCode != '_') {
                fillCode = nextValidCode;
                fillIndex = nextValidIndex;
            } else if (lastValidIndex >= 0 && nextValidIndex == -1 && lastValidCode != '_') {
                fillCode = lastValidCode;
                fillIndex = lastValidIndex;
            }

            if (fillCode != '_') {
                for (int j = lp.getFirst(); j <= lp.getSecond(); j++) {
                    KmerAnnotation ka = smoothedAnnotatedContig.get(j);
                    ka.setCode(fillCode);
                    ka.setBackground(smoothedAnnotatedContig.get(fillIndex).getBackground());
                    ka.setIntervals(selectInterval(annotations.get(j).get(ka.getBackground()), smoothedAnnotatedContig.get(fillIndex).getIntervals()));
                    ka.setSmoothed(true);

                    smoothedAnnotatedContig.set(j, ka);
                }
            }
        }

        return smoothedAnnotatedContig;
    }

    private List<String> splitOnRepeatedCharacters(String annotatedContig, char... splits) {
        Set<Integer> splitPositions = new TreeSet<>();
        for (char split : splits) {
            for (int i = 0; i < annotatedContig.length(); i++) {
                if (annotatedContig.charAt(i) == split) {
                    int j;
                    for (j = i + 1; j < annotatedContig.length() && annotatedContig.charAt(j) == split; j++) {}

                    splitPositions.add(i);
                    splitPositions.add(j);
                }
            }
        }

        List<String> pieces = new ArrayList<>();
        int prevPos = 0;
        for (int nextPos : splitPositions) {
            pieces.add(annotatedContig.substring(prevPos, nextPos));
            prevPos = nextPos;
        }
        pieces.add(annotatedContig.substring(prevPos, annotatedContig.length()));

        return pieces;
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

            //String[] pieces = annotatedContig.split("((?<=\\.+)|(?=\\.+))");
            //annotatedContigs.add(Arrays.asList(pieces));

            List<String> pieces = splitOnRepeatedCharacters(annotatedContig, '.', '_');
            annotatedContigs.add(pieces);

            log.info("  {}", annotatedContig);
            for (int i = 0; i < pieces.size(); i++) {
                log.info("  {} {} {}", background, i, pieces.get(i));
            }
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

        final String alphabet = "0123456789@ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

        IntervalTreeMap<String> itm = new IntervalTreeMap<>();
        for (Map<String, String> m : annotations) {
            int index = Integer.valueOf(m.get("index"));
            if (index == 5596) {
                log.info("H!");
            }

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
