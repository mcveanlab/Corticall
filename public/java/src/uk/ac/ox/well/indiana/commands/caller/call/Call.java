package uk.ac.ox.well.indiana.commands.caller.call;

import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.io.File;
import java.io.PrintStream;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
        List<Integer> recruitColors = GRAPH.getColorsForSampleNames(new ArrayList<>(LOOKUPS.keySet()));

        Map<String, List<Map<String, String>>> contigs = loadAnnotations();

        for (String contigName : contigs.keySet()) {
            callVariants(contigName, contigs.get(contigName));
        }
    }

    private void callVariants(String contigName, List<Map<String, String>> annotations) {
        for (String background : LOOKUPS.keySet()) {
            String annotatedContig = annotateContig(annotations, background);

            if (isNAHR(annotatedContig)) {
                log.info("nahr: {} {}", contigName, annotatedContig);
            }
        }
    }

    private boolean isNAHR(String annotatedContig) {
        final String flankingNovelRegex = "(([^_\\.])\\2+)_*(\\.+)_*(([^_\\.])\\5+)";
        Pattern flankingNovelPattern = Pattern.compile(flankingNovelRegex);

        final String novelRegex = "(\\.+)";
        Pattern novelPattern = Pattern.compile(novelRegex);

        Matcher flankingNovelMatcher = flankingNovelPattern.matcher(annotatedContig);
        int numSwitches = 0;
        if (flankingNovelMatcher.find()) {
            do {
                if (!flankingNovelMatcher.group(2).equals(flankingNovelMatcher.group(5))) {
                    numSwitches++;
                }
            } while (flankingNovelMatcher.find(flankingNovelMatcher.start(3)));
        }

        Matcher novelMatcher = novelPattern.matcher(annotatedContig);
        int numNovelRuns = 0;
        while (novelMatcher.find()) {
            numNovelRuns++;
        }

        return (numSwitches > 0 && numNovelRuns > 1);
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
