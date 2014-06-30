package uk.ac.ox.well.indiana.attic.analyses.nahr;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class LearnContigMetricThresholds extends Module {
    @Argument(fullName="contigMetrics", shortName="cm", doc="Contig metrics")
    public File CM;

    @Argument(fullName="trainingContigs", shortName="tc", doc="Training contigs (BAM)")
    public SAMFileReader TRAINING;

    @Argument(fullName="sampleName", shortName="sn", doc="Sample name")
    public String SAMPLE_NAME;

    @Argument(fullName="metricsBoolean", shortName="mb", doc="Boolean metrics")
    public ArrayList<String> METRICS_BOOLEAN;

    @Argument(fullName="metricsInteger", shortName="mi", doc="Integer metrics")
    public ArrayList<String> METRICS_INTEGER;

    @Output
    public PrintStream out;

    @Output(fullName="contigsOut", shortName="cout", doc="Contigs out")
    public PrintStream cout;

    private List<Map<String, String>> subsetBoolean(List<Map<String, String>> metrics, String metric, int valueDesired) {
        List<Map<String, String>> selectedMetrics = new ArrayList<Map<String, String>>();

        for (Map<String, String> te : metrics) {
            int value = Integer.valueOf(te.get(metric));

            if (valueDesired == value) {
                selectedMetrics.add(te);
            }
        }

        return selectedMetrics;
    }

    private List<Map<String, String>> subsetInteger(List<Map<String, String>> metrics, String metric, int valueDesired) {
        List<Map<String, String>> selectedMetrics = new ArrayList<Map<String, String>>();

        for (Map<String, String> te : metrics) {
            int value = Integer.valueOf(te.get(metric));

            if (valueDesired <= value) {
                selectedMetrics.add(te);
            }
        }

        return selectedMetrics;
    }

    private int numTrainingRecovered(List<Map<String, String>> metrics, Set<String> trainingContigNames) {
        int numRecovered = 0;

        for (Map<String, String> te : metrics) {
            String contigName = te.get("contigName");
            if (trainingContigNames.contains(contigName)) {
                numRecovered++;
            }
        }

        return numRecovered;
    }

    @Override
    public void execute() {
        log.info("Processing training data...");
        Set<String> trainingContigNames = new HashSet<String>();

        for (SAMRecord trainingContig : TRAINING) {
            if (trainingContig.getReadGroup().getSample().equalsIgnoreCase(SAMPLE_NAME)) {
                trainingContigNames.add(trainingContig.getReadName());
            }
        }

        log.info("  found {} training contigs for {}", trainingContigNames.size(), SAMPLE_NAME);

        if (trainingContigNames.size() > 0) {
            log.info("Processing annotated contigs...");
            List<Map<String, String>> metrics = new ArrayList<Map<String, String>>();
            Map<String, Integer> metricMaxValue = new HashMap<String, Integer>();

            TableReader tr = new TableReader(CM);
            for (Map<String, String> te : tr) {
                metrics.add(te);

                for (String metric : te.keySet()) {
                    try {
                        int value = Integer.valueOf(te.get(metric));

                        if (!metricMaxValue.containsKey(metric) || value > metricMaxValue.get(metric)) {
                            metricMaxValue.put(metric, value);
                        }
                    } catch (NumberFormatException e) {}
                }
            }

            log.info("Scanning parameter space for good metric thresholds...");
            List<String> allMetrics = new ArrayList<String>();
            allMetrics.addAll(METRICS_BOOLEAN);
            allMetrics.addAll(METRICS_INTEGER);

            Map<String, Integer> finalValues = new LinkedHashMap<String, Integer>();

            for (String metric : allMetrics) {
                if (METRICS_BOOLEAN.contains(metric)) {
                    List<Map<String, String>> submetrics0 = subsetBoolean(metrics, metric, 0);
                    int numSelected0 = submetrics0.size();
                    int numRecovered0 = numTrainingRecovered(submetrics0, trainingContigNames);

                    List<Map<String, String>> submetrics1 = subsetBoolean(metrics, metric, 1);
                    int numSelected1 = submetrics1.size();
                    int numRecovered1 = numTrainingRecovered(submetrics1, trainingContigNames);

                    finalValues.put(metric, numRecovered0 > numRecovered1 ? 0 : 1);

                    if (numRecovered0 > numRecovered1) {
                        log.info("  {}: value={} selected={} recovered={}", metric, 0, numSelected0, numRecovered0);
                    } else {
                        log.info("  {}: value={} selected={} recovered={}", metric, 1, numSelected1, numRecovered1);
                    }
                } else if (METRICS_INTEGER.contains(metric)) {
                    int finalValue = 0;
                    int numSelected = metrics.size();
                    int numRecovered = numTrainingRecovered(metrics, trainingContigNames);

                    for (int v = 0; v < metricMaxValue.get(metric); v++) {
                        List<Map<String, String>> submetrics = subsetInteger(metrics, metric, v);
                        int newNumSelected = submetrics.size();
                        int newNumRecovered = numTrainingRecovered(submetrics, trainingContigNames);

                        //log.info("  {}: value={} selected={} recovered={}", metric, v, newNumSelected, newNumRecovered);

                        if (newNumRecovered == numRecovered) {
                            finalValue = v;
                            numSelected = newNumSelected;
                            numRecovered = newNumRecovered;
                        } else {
                            break;
                        }
                    }

                    finalValues.put(metric, finalValue == 0 && numSelected == metrics.size() ? 0 : finalValue - 1);

                    log.info("  {}: value={} selected={} recovered={}", metric, finalValues.get(metric), numSelected, numRecovered);
                }
            }

            List<Map<String, String>> subset = metrics;

            out.println("sample\tmetric\tvalue");
            for (String metric : finalValues.keySet()) {
                out.println(SAMPLE_NAME + "\t" + metric + "\t" + finalValues.get(metric));

                if (METRICS_BOOLEAN.contains(metric)) {
                    subset = subsetBoolean(subset, metric, finalValues.get(metric));
                } else {
                    subset = subsetInteger(subset, metric, finalValues.get(metric));
                }
            }

            for (Map<String, String> te : subset) {
                cout.println(te.get("contigName"));
            }

            log.info("Final: selected={} recovered={}", subset.size(), numTrainingRecovered(subset, trainingContigNames));
        }
    }
}
