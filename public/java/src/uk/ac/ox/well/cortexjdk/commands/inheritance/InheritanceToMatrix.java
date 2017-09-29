package uk.ac.ox.well.cortexjdk.commands.inheritance;

import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableWriter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by kiran on 20/08/2017.
 */
public class InheritanceToMatrix extends Module {
    @Argument(fullName="track", shortName="t", doc="Track")
    public File TRACK;

    @Argument(fullName="parent0", shortName="p0", doc="Parent 0")
    public String PARENT0;

    @Argument(fullName="parent1", shortName="p1", doc="Parent 1")
    public String PARENT1;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, Map<String, List<Pair<Integer, Boolean>>>> tracks = loadTracks(TRACK);
        Map<String, Map<String, List<Pair<Integer, Boolean>>>> smoothedTracks = smoothTracks(tracks);

        printMatrix(smoothedTracks);
    }

    private void printMatrix(Map<String, Map<String, List<Pair<Integer, Boolean>>>> simplifiedTracks) {
        log.info("Printing tracks...");

        TableWriter tw = new TableWriter(out);

        String firstSampleName = simplifiedTracks.keySet().iterator().next();

        for (String chrom : simplifiedTracks.get(firstSampleName).keySet()) {
            int length = simplifiedTracks.get(firstSampleName).get(chrom).size();

            for (int i = 0; i < length; i++) {
                Map<String, String> te = new LinkedHashMap<>();

                for (String sampleName : simplifiedTracks.keySet()) {
                    te.put(sampleName, simplifiedTracks.get(sampleName).get(chrom).get(0).getSecond() ? "1" : "0");
                }

                tw.addEntry(te);
            }
        }
    }

    @NotNull
    private Map<String, Map<String, List<Pair<Integer, Boolean>>>> loadTracks(File trackFile) {
        TableReader tr = new TableReader(trackFile);

        Map<String, Map<String, List<Pair<Integer, Boolean>>>> tracks = new HashMap<>();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing track records")
                .message("records processed")
                .maxRecord(tr.size())
                .make(log);

        for (Map<String, String> te : tr) {
            String chrom = te.get("chrom");
            int pos = Integer.valueOf(te.get("pos"));
            String type = te.get("type");
            //int covParent = Integer.valueOf(te.get("cov_parent"));
            //String alleles = te.get("alleles");

            if (type.equals("SNP")) {
                for (String sampleName : te.keySet()) {
                    if (!sampleName.equals("chrom") && !sampleName.equals("pos") && !sampleName.equals("type") && !sampleName.equals("cov_parent") && !sampleName.contains("alleles")) {
                        String[] pieces = te.get(sampleName).split(":");
                        boolean parentCode = pieces[0].equals(PARENT1);

                        Pair<Integer, Boolean> p = new Pair<>(pos, parentCode);

                        if (!tracks.containsKey(sampleName)) { tracks.put(sampleName, new TreeMap<>()); }
                        if (!tracks.get(sampleName).containsKey(chrom)) { tracks.get(sampleName).put(chrom, new ArrayList<>()); }

                        tracks.get(sampleName).get(chrom).add(p);
                    }
                }
            }

            pm.update();
        }
        return tracks;
    }

    Map<String, Map<String, List<Pair<Integer, Boolean>>>> simplifyTracks(Map<String, Map<String, List<Pair<Integer, Boolean>>>> tracks) {
        log.info("Simplifying tracks...");

        Map<String, Map<String, List<Pair<Integer, Boolean>>>> simplifiedTracks = new HashMap<>();

        for (String sampleName : tracks.keySet()) {
            simplifiedTracks.put(sampleName, new TreeMap<>());

            for (String chrom : tracks.get(sampleName).keySet()) {
                List<Pair<Integer, Boolean>> l = tracks.get(sampleName).get(chrom);
                List<Pair<Integer, Boolean>> t = new ArrayList<>();

                t.add(l.get(0));

                for (int i = 1; i < l.size(); i++) {
                    Pair<Integer, Boolean> p = l.get(i);

                    if (t.get(t.size() - 1).getSecond() != p.getSecond() || i == l.size() - 1) {
                        t.add(p);
                    }
                }

                simplifiedTracks.get(sampleName).put(chrom, t);
            }
        }

        return simplifiedTracks;
    }

    Map<String, Map<String, List<Pair<Integer, Boolean>>>> smoothTracks(Map<String, Map<String, List<Pair<Integer, Boolean>>>> tracks) {
        log.info("Smoothed tracks...");

        Map<String, Map<String, List<Pair<Integer, Boolean>>>> smoothTracks = new HashMap<>();

        for (String sampleName : tracks.keySet()) {
            smoothTracks.put(sampleName, new TreeMap<>());

            for (String chrom : tracks.get(sampleName).keySet()) {
                List<Pair<Integer, Boolean>> l = tracks.get(sampleName).get(chrom);
                List<Pair<Integer, Boolean>> t = new ArrayList<>();

                t.add(l.get(0));

                for (int i = 1; i < l.size() - 1; i++) {
                    Pair<Integer, Boolean> pm1 = l.get(i-1);
                    Pair<Integer, Boolean> pi = l.get(i);
                    Pair<Integer, Boolean> pp1 = l.get(i+1);

                    if (pm1.getSecond() == pp1.getSecond() && pm1.getSecond() != pi.getSecond()) {
                        t.add(new Pair<>(pi.getFirst(), pm1.getSecond()));
                    } else {
                        t.add(pi);
                    }
                }

                t.add(l.get(l.size() - 1));

                smoothTracks.get(sampleName).put(chrom, t);
            }
        }

        return smoothTracks;
    }

    private Map<String, PrintStream> initializeOutputs(File prefix, Set<String> headers) {
        Map<String, PrintStream> outs = new HashMap<>();

        for (String header : headers) {
            if (!header.equals("chrom") && !header.equals("pos") && !header.equals("type") && !header.equals("cov_parent") && !header.contains("alleles")) {
                File outFile = new File(prefix.getAbsolutePath() + "." + header + ".circos.txt");
                try {
                    PrintStream out = new PrintStream(outFile);

                    outs.put(header, out);
                } catch (FileNotFoundException e) {
                    throw new CortexJDKException("Could not open file '" + outFile.getAbsolutePath() + "' for writing.");
                }
            }
        }

        return outs;
    }
}
