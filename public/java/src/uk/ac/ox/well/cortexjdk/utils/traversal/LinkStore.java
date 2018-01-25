package uk.ac.ox.well.cortexjdk.utils.traversal;

import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexJunctionsRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinksRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.util.*;

/**
 * Created by kiran on 24/07/2017.
 */
public class LinkStore {
    private Map<String, List<LinkStoreElement>> linkElements = new HashMap<>();

    public void add(CortexByteKmer curKmer, CortexLinksRecord clr, boolean goForward, String linkSource) {
        boolean recordOrientationMatchesKmer = clr.getKmerAsByteKmer().equals(curKmer);

        List<CortexJunctionsRecord> junctions = new ArrayList<>();
        junctions.addAll(clr.getJunctions());

        for (CortexJunctionsRecord cjr : junctions) {
            //boolean linkGoesForward = recordOrientationMatchesKmer ? cjr.isForward() : !cjr.isForward();
            boolean linkGoesForward = recordOrientationMatchesKmer == cjr.isForward();
            String junctionList = linkGoesForward ? cjr.getJunctions() : SequenceUtils.complement(cjr.getJunctions());

            if (linkGoesForward == goForward) {
                if (!linkElements.containsKey(junctionList)) {
                    linkElements.put(junctionList, new ArrayList<>());
                }

                linkElements.get(junctionList).add(new LinkStoreElement(junctionList, 0, 0, linkSource));
            }
        }
    }

    public void incrementAges() {
        for (String junctionList : linkElements.keySet()) {
            for (LinkStoreElement lse : linkElements.get(junctionList)) {
                lse.incrementAge();
            }
        }
    }

    public void imcrementPositions() {
        for (String junctionList : linkElements.keySet()) {
            for (LinkStoreElement lse : linkElements.get(junctionList)) {
                lse.incrementPos();
            }
        }
    }

    public int numNewPaths() {
        int numNewPaths = 0;
        for (String junctionList : linkElements.keySet()) {
            for (LinkStoreElement lse : linkElements.get(junctionList)) {
                if (lse.getAge() == 0) {
                    numNewPaths++;
                }
            }
        }

        return numNewPaths;
    }

    /*
    private void expireRecords() {
        Set<LinkStoreElement> toRemove = new HashSet<>();

        List<LinkStoreElement> lses = new ArrayList<>();
        for (String junctionList : linkElements.keySet()) {
            for (LinkStoreElement lse : linkElements.get(junctionList)) {
                lses.add(lse);
            }
        }

        lses.sort((o1, o2) -> {
            if (o1.getAge() != o2.getAge()) { return o1.getAge() > o2.getAge() ? -1 : 1; }
            if (o1.length() == o2.length()) { return o2.getJunctionList().compareTo(o1.getJunctionList()); }
            return o1.length() < o2.length() ? -1 : 1;
        });

        for (LinkStoreElement lse : lses) {
            int pos = lse.getPos();
            String junctionList = lse.getJunctionList();

            if (lse.getPos() + 1 >= junctionList.length() || !junctionList.substring(pos, pos + 1).equals(choice)) {
                toRemove.add(lse);
            } else if (junctionList.substring(pos, pos + 1).equals(choice)) {
                lse.incrementPos();
            }
        }
    }
    */

    private void incrementPositionsAndExpire(String choice) {
        Set<LinkStoreElement> toRemove = new HashSet<>();

        List<LinkStoreElement> lses = new ArrayList<>();
        for (String junctionList : linkElements.keySet()) {
            for (LinkStoreElement lse : linkElements.get(junctionList)) {
                lses.add(lse);
            }
        }

        lses.sort((o1, o2) -> {
            if (o1.getAge() != o2.getAge()) { return o1.getAge() > o2.getAge() ? -1 : 1; }
            if (o1.length() == o2.length()) { return o2.getJunctionList().compareTo(o1.getJunctionList()); }
            return o1.length() < o2.length() ? -1 : 1;
        });

        for (LinkStoreElement lse : lses) {
            int pos = lse.getPos();
            String junctionList = lse.getJunctionList();

            if (lse.getPos() + 1 >= junctionList.length() || !junctionList.substring(pos, pos + 1).equals(choice)) {
                toRemove.add(lse);
            } else if (junctionList.substring(pos, pos + 1).equals(choice)) {
                lse.incrementPos();
            }
        }

        for (LinkStoreElement lse : toRemove) {
            linkElements.get(lse.getJunctionList()).remove(lse);

            if (linkElements.get(lse.getJunctionList()).size() == 0) {
                linkElements.remove(lse.getJunctionList());
            }
        }
    }

    private String getOldestLink() {
        int age = Integer.MIN_VALUE;

        for (String junctionList : linkElements.keySet()) {
            for (LinkStoreElement lse : linkElements.get(junctionList)) {
                if (lse.getAge() > age) {
                    age = lse.getAge();
                }
            }
        }

        Set<LinkStoreElement> oldestLinks = new LinkedHashSet<>();
        for (String junctionList : linkElements.keySet()) {
            for (LinkStoreElement lse : linkElements.get(junctionList)) {
                if (lse.getAge() == age) {
                    oldestLinks.add(lse);
                }
            }
        }

        Set<String> junctionChoices = new LinkedHashSet<>();
        for (LinkStoreElement lse : oldestLinks) {
            if (lse.getPos() + 1 <= lse.length()) {
                junctionChoices.add(lse.getJunctionList().substring(lse.getPos(), lse.getPos() + 1));
            }
        }

        return junctionChoices.size() == 1 ? oldestLinks.iterator().next().getJunctionList() : null;
    }

    public Pair<String, Set<String>> getNextJunctionChoice() {
        String junctionList = getOldestLink();
        String choice = null;
        Set<String> junctionSources = new TreeSet<>();

        if (junctionList != null) {
            for (LinkStoreElement lse : linkElements.get(junctionList)) {
                int pos = lse.getPos();

                choice = junctionList.substring(pos, pos + 1);

                for (String jl : linkElements.keySet()) {
                    if (pos < jl.length() && jl.substring(pos, pos + 1).equals(choice)) {
                        junctionSources.add(lse.getSource());
                    }
                }
            }

            incrementPositionsAndExpire(choice);
        }

        return new Pair<>(choice, junctionSources);
    }

    public boolean isActive() {
        return linkElements.size() > 0;
    }

    public int size() {
        int numElements = 0;

        for (String junctionList : linkElements.keySet()) {
            numElements += linkElements.get(junctionList).size();
        }

        return numElements;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();

        List<LinkStoreElement> junctionLists = new ArrayList<>();

        for (String jl : linkElements.keySet()) {
            junctionLists.addAll(linkElements.get(jl));
        }

        junctionLists.sort((o1, o2) -> {
            if (o1.getAge() != o2.getAge()) {
                return o1.getAge() > o2.getAge() ? -1 : 1;
            }

            if (o1.length() == o2.length()) { return o2.getJunctionList().compareTo(o1.getJunctionList()); }
            return o1.length() < o2.length() ? -1 : 1;
        });

        int numCurr = junctionLists.size();

        sb.append("\n");
        sb.append("                             num_curr: ").append(numCurr);

        for (LinkStoreElement lse : junctionLists) {
            if (sb.length() > 0) {
                sb.append("\n");
            }

            sb.append("                              ")
                    .append(SequenceUtils.complement(lse.getJunctionList()))
                    .append(" [")
                    .append(lse.getPos())
                    .append("/")
                    .append(lse.length())
                    .append("] age: ")
                    .append(lse.getAge());
        }

        sb.append("\n");
        sb.append("                             num_counter: 0\n");

        return sb.toString();
    }
}
