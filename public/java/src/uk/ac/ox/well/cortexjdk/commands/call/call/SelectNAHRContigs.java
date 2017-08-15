package uk.ac.ox.well.cortexjdk.commands.call.call;

import com.google.common.base.Joiner;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by kiran on 09/08/2017.
 */
public class SelectNAHRContigs extends Module {
    @Argument(fullName="annotations", shortName="a", doc="Annotated contigs")
    public File ANNOTATIONS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, List<Map<String, String>>> allAnnotations = loadAnnotations();

        int nahrEventId = 0;

        for (String contigName : allAnnotations.keySet()) {
            StringBuilder sb = new StringBuilder();

            for (Map<String, String> m : allAnnotations.get(contigName)) {
                String code = m.get("code");

                sb.append(code);
            }

            if (isNahrEvent(sb.toString())) {
                log.info("nahr: {} {}", contigName, sb);

                for (int i = 0; i < sb.length(); i++) {
                    Map<String, String> m = allAnnotations.get(contigName).get(i);

                    if (m.get("code").equals(".")) {
                        String prevInterval = null;
                        String nextInterval = null;

                        for (int j = i - 1; j >= 0; j--) {
                            Map<String, String> mj = allAnnotations.get(contigName).get(j);
                            if (!mj.get("intervals").equals("NA")) {
                                prevInterval = mj.get("intervals");
                                break;
                            }
                        }

                        for (int j = i + 1; j < sb.length(); j++) {
                            Map<String, String> mj = allAnnotations.get(contigName).get(j);
                            if (!mj.get("intervals").equals("NA")) {
                                nextInterval = mj.get("intervals");
                                i = j;
                                break;
                            }
                        }

                        if (prevInterval != null && nextInterval != null) {
                            if (prevInterval.contains("SI_SI") || prevInterval.contains("_NIH_CSHL")) {
                                prevInterval = prevInterval.replaceAll("HB3|DD2|7G8|GB4|803", "3D7");
                                prevInterval = prevInterval.replaceAll("_SI_SI|_NIH_CSHL", "");
                                prevInterval = prevInterval.replaceFirst(":", "_v3:");
                            }

                            if (nextInterval.contains("SI_SI") || nextInterval.contains("_NIH_CSHL")) {
                                nextInterval = nextInterval.replaceAll("HB3|DD2|7G8|GB4|803", "3D7");
                                nextInterval = nextInterval.replaceAll("_SI_SI|_NIH_CSHL", "");
                                nextInterval = nextInterval.replaceFirst(":", "_v3:");
                            }

                            String[] pi = prevInterval.split("[:-]");
                            String[] ni = nextInterval.split("[:-]");

                            out.println(Joiner.on(" ").join(pi[0], pi[1], pi[2], ni[0], ni[1], ni[2], "thickness=" + nahrEventId));

                            //log.info("  {} {}", prevInterval, nextInterval);
                        }
                    }
                }

                nahrEventId++;
            }
        }
    }

    private boolean isNahrEvent(String annotatedContig) {
        int numTemplateSwitches = numTemplateSwitches(annotatedContig);
        int numNovelRuns = numNovelRuns(annotatedContig);

        if (numTemplateSwitches >= 2 && numNovelRuns >= 1) {
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

    private Map<String, List<Map<String, String>>> loadAnnotations() {
        TableReader tr = new TableReader(ANNOTATIONS);

        Map<String, List<Map<String, String>>> contigs = new TreeMap<>();

        for (Map<String, String> m : tr) {
            if (!contigs.containsKey(m.get("contigName"))) {
                contigs.put(m.get("contigName"), new ArrayList<>());
            }
            contigs.get(m.get("contigName")).add(m);
        }

        return contigs;
    }
}
