package uk.ac.ox.well.cortexjdk.commands.discover.call;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.util.*;

public class IntervalCombiner {
    private IntervalCombiner() {}

    public static List<Pair<String, Interval>> getIntervals(List<CortexVertex> ws, Map<String, IndexedReference> refs, int locationWindow, int maxStatesPerBackground) {
        Map<String, IntervalTreeMap<Interval>> itm = new HashMap<>();

        for (String name : refs.keySet()) {
            itm.put(name, new IntervalTreeMap<>());

            for (int i = 0; i < ws.size(); i++) {
                CortexVertex cv = ws.get(i);

                for (Interval it : refs.get(name).find(cv.getKmerAsString())) {
                    Interval pit = new Interval(it.getContig(), it.getStart() - locationWindow, it.getEnd() + locationWindow, it.isNegativeStrand(), it.getName());

                    if (itm.get(name).containsOverlapping(pit) || itm.get(name).containsContained(pit)) {
                        int newStart = pit.getStart();
                        int newEnd = pit.getEnd();

                        Set<Interval> toRemove = new HashSet<>();
                        for (Collection<Interval> intervals : Arrays.asList(itm.get(name).getOverlapping(pit), itm.get(name).getContained(pit))) {
                            for (Interval oit : intervals) {
                                if (pit.isNegativeStrand() == oit.isNegativeStrand()) {
                                    if (oit.getStart() < newStart) {
                                        newStart = oit.getStart();
                                    }
                                    if (oit.getEnd() > newEnd) {
                                        newEnd = oit.getEnd();
                                    }

                                    toRemove.add(oit);
                                }
                            }
                        }

                        pit = new Interval(it.getContig(), newStart, newEnd, it.isNegativeStrand(), it.getName());

                        for (Interval rit : toRemove) {
                            itm.get(name).remove(rit);
                        }
                    }

                    itm.get(name).put(pit, pit);
                }
            }
        }

        List<Pair<String, Interval>> allStates = new ArrayList<>();
        for (String s : itm.keySet()) {
            Set<Interval> sortedStates = new TreeSet<>(Comparator.comparingInt(Interval::length).reversed());

            for (Interval it : itm.get(s).values()) {
                if (it.getStart() < 0) {
                    it = new Interval(it.getContig(), 1, it.getEnd(), it.isNegativeStrand(), it.getName());
                }

                int maxLength = refs.get(s).getReferenceSequence().getSequenceDictionary().getSequence(it.getContig()).getSequenceLength();
                if (it.getEnd() >= maxLength) {
                    it = new Interval(it.getContig(), it.getStart(), maxLength - 1, it.isNegativeStrand(), it.getName());
                }

                //String seq = refs.get(s).find(it);
                //String targetName = s + ":" + it.getContig() + ":" + it.getStart() + "-" + it.getEnd() + ":" + (it.isPositiveStrand() ? "+" : "-");
                //targets.put(targetName, seq);

                sortedStates.add(it);
            }

            Iterator<Interval> it = sortedStates.iterator();
            for (int i = 0; i < maxStatesPerBackground && it.hasNext(); i++) {
                allStates.add(new Pair<>(s, it.next()));
            }
        }

        return allStates;
    }
}
