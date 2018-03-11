package uk.ac.ox.well.cortexjdk.commands.discover.call;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.util.*;

import static uk.ac.ox.well.cortexjdk.commands.discover.call.BackgroundInference.BackgroundGroup.CHILD;
import static uk.ac.ox.well.cortexjdk.commands.discover.call.BackgroundInference.BackgroundGroup.PARENT;
import static uk.ac.ox.well.cortexjdk.commands.discover.call.BackgroundInference.BackgroundGroup.REF;

public class BackgroundInference {
    public enum BackgroundGroup { CHILD, PARENT, REF }

    private List<Integer> indices = new ArrayList<>();
    private List<CortexVertex> vertices = new ArrayList<>();
    private List<String> names = new ArrayList<>();
    private List<Set<Interval>> locations = new ArrayList<>();
    private List<Boolean> presence = new ArrayList<>();
    private List<Boolean> novelty = new ArrayList<>();
    private List<BackgroundGroup> groups = new ArrayList<>();

    private int locationWindow = 200;
    private Map<BackgroundGroup, float[][][]> transitionMatrices = null;

    public BackgroundInference() {}

    public BackgroundInference(int locationWindow) {
        this.locationWindow = locationWindow;
    }

    public void add(int index, CortexVertex vertex, String name, Set<Interval> locations, boolean isPresent, boolean isNovel, BackgroundGroup group) {
        this.indices.add(index);
        this.vertices.add(vertex);
        this.names.add(name);
        this.locations.add(locations);
        this.presence.add(isPresent);
        this.novelty.add(isNovel);
        this.groups.add(group);
    }

    public List<Pair<String, Interval>> getStates(BackgroundGroup group, int maxStatesPerBackground) {
        Map<String, IntervalTreeMap<Interval>> itm = new HashMap<>();

        for (int i = 0; i < indices.size(); i++) {
            if (groups.get(i) == group) {
                String name = names.get(i);

                if (!itm.containsKey(name)) {
                    itm.put(name, new IntervalTreeMap<>());
                }

                for (Interval it : locations.get(i)) {
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
            //for (Interval it : itm.get(s).values()) {
                //allStates.add(new Pair<>(s, it));
                //sortedStates.add(it);
            //}

            sortedStates.addAll(itm.get(s).values());

            Iterator<Interval> it = sortedStates.iterator();
            for (int i = 0; i < maxStatesPerBackground && it.hasNext(); i++) {
                allStates.add(new Pair<>(s, it.next()));
            }
        }

        return allStates;
    }

    public Triple<List<String>, double[][][], List<String>> prepare(BackgroundGroup group) {
        int indexMin = Integer.MAX_VALUE;
        int indexMax = Integer.MIN_VALUE;
        Map<String, IntervalTreeMap<Interval>> itm = new HashMap<>();

        for (int i = 0; i < indices.size(); i++) {
            if (groups.get(i) == group) {
                if (indices.get(i) < indexMin) { indexMin = indices.get(i); }
                if (indices.get(i) > indexMax) { indexMax = indices.get(i); }

                String name = names.get(i);

                if (!itm.containsKey(name)) {
                    itm.put(name, new IntervalTreeMap<>());
                }

                for (Interval it : locations.get(i)) {
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

        boolean[] isNovel = new boolean[indexMax + 1];
        for (int i = 0; i < indices.size(); i++) {
            if (groups.get(i) == CHILD) {
                isNovel[indices.get(i)] = novelty.get(i);
            }
        }

        Set<String> stateSet = new LinkedHashSet<>();

        for (String s : itm.keySet()) {
            Set<Interval> sortedIntervals = new TreeSet<>(Comparator.comparingInt(Interval::length).reversed());

            for (Interval it : new TreeSet<>(itm.get(s).keySet())) {
                Interval pit = itm.get(s).get(it);

                sortedIntervals.add(pit);
            }

            Iterator<Interval> iter = sortedIntervals.iterator();
            for (int i = 0; i < Math.min(sortedIntervals.size(), 50); i++) {
                Interval pit = iter.next();

                String state = String.format("%s:%s:%d-%d:%s", s, pit.getContig(), pit.getStart(), pit.getEnd(), pit.isPositiveStrand() ? "+" : "-");
                stateSet.add(state);
            }

            stateSet.add(s);
        }

        List<String> stateList = new ArrayList<>(stateSet);

        double[][][] a = new double[indexMax - indexMin + 1][stateList.size()][stateList.size()];

        for (int index = 0; index < a.length; index++) {
            for (int si1 = 0; si1 < stateList.size(); si1++) {
                for (int si2 = 0; si2 < stateList.size(); si2++) {
                    a[index][si1][si2] = 1.0 / stateList.size();
                }
            }
        }

        for (int i = 0; i < indices.size(); i++) {
            if (groups.get(i) == group) {
                String s = names.get(i);
                int index = indices.get(i) - indexMin;
                boolean isPresent = presence.get(i);

                Set<String> states = new HashSet<>();

                Set<Interval> intervals = locations.get(i);
                Set<Interval> refinedIntervals = new HashSet<>();

                if (intervals.size() > 0) {
                    for (Interval it : intervals) {
                        Interval rit = new Interval(it.getContig(), it.getStart(), it.getEnd());
                        Collection<Interval> pits = itm.get(s).getOverlapping(rit);

                        refinedIntervals.addAll(pits);
                    }
                }

                if (refinedIntervals.size() > 0) {
                    for (Interval rit : refinedIntervals) {
                        String state = String.format("%s:%s:%d-%d:%s", s, rit.getContig(), rit.getStart(), rit.getEnd(), rit.isPositiveStrand() ? "+" : "-");

                        states.add(state);
                    }
                } else {
                    states.add(s);
                }

                for (String state : states) {
                    int stateIndex = stateList.indexOf(state);

                    if (stateIndex >= 0) {
                        double stayTheSameState = 100000;
                        double switchMyParentBase = 1;
                        double switchMyParentState = 1;
                        double switchOtherParentState = 1;

                        if (!isPresent && index > 0) {
                            if ((!isNovel[index - 1] && isNovel[index]) || (isNovel[index - 1] && !isNovel[index])) {
                                stayTheSameState = 1;
                                switchMyParentBase = 10;
                                switchMyParentState = 10;
                                switchOtherParentState = 10;
                            } else if (isNovel[index - 1] && isNovel[index]) {
                                stayTheSameState = 1000;
                                switchMyParentBase = 1;
                                switchMyParentState = 1;
                                switchOtherParentState = 1;
                            } else if (!isNovel[index - 1] && !isNovel[index]) {
                                stayTheSameState = 1;
                                switchMyParentBase = 100;
                                switchMyParentState = 1;
                                switchOtherParentState = 1;
                            }
                        }

                        double norm = 0;
                        for (int si2 = 0; si2 < stateList.size(); si2++) {
                            if (stateIndex == si2) { a[index][stateIndex][si2] = stayTheSameState; }
                            else {
                                if (stateList.get(si2).equals(s)) {
                                    a[index][stateIndex][si2] = switchMyParentBase;
                                } else if (stateList.get(si2).startsWith(s)) {
                                    a[index][stateIndex][si2] = switchMyParentState;
                                } else {
                                    a[index][stateIndex][si2] = switchOtherParentState;
                                }
                            }

                            norm += a[index][stateIndex][si2];
                        }

                        for (int si2 = 0; si2 < stateList.size(); si2++) {
                            a[index][stateIndex][si2] /= norm;
                        }

                        //a[index][stateIndex][stateIndex] = isPresent ? 0.99999f : 2.0 * (1.0 / stateList.size());

                        /*
                        a[index][stateIndex][stateIndex] = isPresent ? 0.99999 : 2.0 * (1.0 / stateList.size());
                        double rem = a[index][stateIndex][stateIndex] / (stateList.size());

                        for (int si1 = 0; si1 < stateList.size(); si1++) {
                            for (int si2 = 0; si2 < stateList.size(); si2++) {
                                if (si1 != si2) {
                                    a[index][si1][si2] = rem;
                                }
                            }
                        }

                        if (index > 0 && isNovel[index] && !isNovel[index - 1]) {
                            a[index][stateIndex][stateIndex] = 0;
                            a[index][stateIndex][stateList.indexOf(s)] = 0.5;

                            rem = a[index][stateIndex][stateList.indexOf(s)] / (stateList.size());

                            for (int si2 = 0; si2 < stateList.size(); si2++) {
                                a[index][stateIndex][si2] = rem;
                            }
                        }
                        */
                    }
                }
            }
        }

        List<String> stateSequence = computeStateSequence(stateList, a);

        return Triple.of(stateList, a, stateSequence);
    }

    public List<String> computeStateSequence(List<String> stateList, double[][][] a) {
        double[][] v = new double[a.length][stateList.size()];
        int[][] t = new int[a.length][stateList.size()];

        for (int k = 0; k < v[0].length; k++) {
            v[0][k] = Math.log(1.0 / ((double) stateList.size()));
            t[0][k] = 0;
        }

        for (int i = 1; i < a.length; i++) {
            for (int k = 0; k < v[0].length; k++) {
                v[i][k] = Math.log10(Double.MIN_VALUE);
            }
        }

        for (int i = 1; i < a.length; i++) {
            for (int j = 0; j < v[0].length; j++) {
                for (int k = 0; k < v[0].length; k++) {
                    double aijk = Math.log10(a[i-1][j][k] + Double.MIN_VALUE);
                    if (v[i][j] < (v[i-1][k]) + aijk) {
                        v[i][j] = (v[i-1][k]) + aijk;
                        t[i][j] = k;
                    }
                }
            }
        }

        /*
        for (int i = 0; i < a.length; i++) {
            List<String> values = new ArrayList<>();

            Arrays.stream(v[i]).forEach(q -> values.add(String.format("%.2f", q)));

            Main.getLogger().info("{} {}", i, Joiner.on(" ").join(values));
        }
        */

        int[] z = new int[a.length];
        String[] x = new String[a.length];

        double q = 0.0;
        for (int k = 0; k < v[0].length; k++) {
            if (v[a.length - 1][k] > q) {
                q = v[a.length - 1][k];
                z[a.length - 1] = k;
            }
        }

        x[a.length - 1] = stateList.get(z[z.length - 1]);

        for (int i = a.length - 1; i >= 1; i--) {
            z[i-1] = t[i][z[i]];
            x[i-1] = stateList.get(z[i-1]);
        }

        return Arrays.asList(x);
    }

    public void weight(List<CortexVertex> vs) {
        if (transitionMatrices == null) {
            for (BackgroundGroup bg : Arrays.asList(PARENT, REF)) {
                //transitionMatrices.put(bg, prepare(bg));
            }
        }
    }
}
