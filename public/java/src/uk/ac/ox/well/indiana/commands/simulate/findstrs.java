package uk.ac.ox.well.indiana.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.PrintStream;
import java.util.*;

public class findstrs extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference FASTA")
    public FastaSequenceFile REFERENCE;

    @Argument(fullName="repeatThreshold", shortName="t", doc="Repeat threshold")
    public Integer REPEAT_THRESHOLD = 3;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableWriter tw = new TableWriter(out);

        log.info("Finding STRs in reference genome...");

        ReferenceSequence rseq;
        while ((rseq = REFERENCE.nextSequence()) != null) {
            String seq = new String(rseq.getBases());
            String[] name = rseq.getName().split("\\s+");

            for (int strLength = 2; strLength <= 5; strLength++) {
                Map<String, IntervalTreeMap<Integer>> strs = new HashMap<String, IntervalTreeMap<Integer>>();

                // Make a list of all kmers in the chromosome and where they come from
                for (int i = 0; i <= seq.length() - strLength; i++) {
                    String kmer = seq.substring(i, i + strLength);

                    if (!strs.containsKey(kmer)) {
                        strs.put(kmer, new IntervalTreeMap<Integer>());
                    }

                    Interval interval = new Interval(name[0], i, i + strLength);
                    strs.get(kmer).put(interval, null);
                }

                // Make a list of kmers that occur less than 2 times in the chromosome
                Set<String> kmersToRemove = new HashSet<String>();
                for (String kmer : strs.keySet()) {
                    if (strs.get(kmer).size() < 2) {
                        kmersToRemove.add(kmer);
                    }
                }

                // Remove aforementioned kmers
                for (String kmer : kmersToRemove) {
                    strs.remove(kmer);
                }

                // Merge adjacent intervals
                for (String kmer : strs.keySet()) {
                    Set<Interval> intervalsToRemove = new HashSet<Interval>();
                    IntervalTreeMap<Integer> expandedMap = new IntervalTreeMap<Integer>();

                    for (Interval interval : strs.get(kmer).keySet()) {
                        Interval expandedInterval = new Interval(interval.getSequence(), interval.getStart(), interval.getEnd());
                        int numRepeats = 1;
                        boolean keepGoing = true;

                        if (!intervalsToRemove.contains(interval)) {
                            do {
                                Interval nextInterval = new Interval(expandedInterval.getSequence(),
                                        expandedInterval.getEnd(),
                                        expandedInterval.getEnd() + strLength);

                                if (strs.get(kmer).containsKey(nextInterval)) {
                                    expandedInterval = new Interval(expandedInterval.getSequence(),
                                            expandedInterval.getStart(),
                                            expandedInterval.getEnd() + strLength);

                                    intervalsToRemove.add(nextInterval);

                                    numRepeats++;
                                } else {
                                    keepGoing = false;
                                }
                            } while (keepGoing);

                            intervalsToRemove.add(interval);
                        }

                        expandedMap.put(expandedInterval, numRepeats);
                    }

                    strs.put(kmer, expandedMap);
                }

                // Retain STRs longer than one repeat
                Map<String, IntervalTreeMap<Integer>> goodStrs = new HashMap<String, IntervalTreeMap<Integer>>();

                for (String kmer : strs.keySet()) {
                    Set<Interval> intervals = strs.get(kmer).keySet();

                    for (Interval interval : intervals) {
                        Set<Integer> usedRepeats = new HashSet<Integer>();
                        for (Integer numRepeats : strs.get(kmer).getContained(interval)) {
                            if (numRepeats >= REPEAT_THRESHOLD && !usedRepeats.contains(numRepeats)) {
                                if (!goodStrs.containsKey(kmer)) {
                                    goodStrs.put(kmer, new IntervalTreeMap<Integer>());
                                }

                                goodStrs.get(kmer).put(interval, numRepeats);

                                usedRepeats.add(numRepeats);
                            }
                        }
                    }
                }

                // Print
                Map<String, Integer> strCounts = new HashMap<String, Integer>();
                for (String kmer : goodStrs.keySet()) {
                    Set<Interval> intervals = goodStrs.get(kmer).keySet();

                    strCounts.put(kmer, intervals.size());

                    for (Interval interval : intervals) {
                        //log.info("    {} {}: {}", kmer, interval, Joiner.on(", ").join(goodStrs.get(kmer).getContained(interval)));

                        int numRepeats = goodStrs.get(kmer).getContained(interval).iterator().next();

                        Map<String, String> te = new LinkedHashMap<String, String>();
                        te.put("strLength", String.valueOf(strLength));
                        te.put("str", kmer);
                        te.put("numRepeats", String.valueOf(numRepeats));
                        te.put("chr", interval.getSequence());
                        te.put("start", String.valueOf(interval.getStart()));
                        te.put("end", String.valueOf(interval.getEnd()));

                        tw.addEntry(te);
                    }
                }

                log.info("  {}: {} {} ({})", name[0], strLength, goodStrs.size(), Joiner.on(", ").withKeyValueSeparator("=").join(strCounts));
            }
        }
    }
}
