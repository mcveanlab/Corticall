package uk.ac.ox.well.indiana.utils.assembly;

import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.util.*;

public class CortexGraphWalker {
    Map<String, CortexRecord> records;

    public CortexGraphWalker(Map<String, CortexRecord> records) {
        this.records = records;
    }

    private List<String> getNextInKmers(String kmer, String edges) {
        List<String> inRawKmers = new ArrayList<String>();

        for (int i = 0; i < 4; i++) {
            if (edges.charAt(i) != '.') {
                String inRawKmer = (edges.charAt(i) + kmer.substring(0, kmer.length() - 1)).toUpperCase();

                inRawKmers.add(inRawKmer);
            }
        }

        return inRawKmers;
    }

    private List<String> getNextInKmers(int color, String kmer) {
        CortexRecord record = records.get(kmer);
        String edges = record.getEdges(color);

        return getNextInKmers(kmer, edges);
    }

    private List<String> getNextOutKmers(String kmer, String edges) {
        List<String> outRawKmers = new ArrayList<String>();

        for (int i = 4; i < 8; i++) {
            if (edges.charAt(i) != '.') {
                String outRawKmer = (kmer.substring(1, kmer.length()) + edges.charAt(i)).toUpperCase();

                outRawKmers.add(outRawKmer);
            }
        }

        return outRawKmers;
    }

    private List<String> getNextOutKmers(int color, String kmer) {
        CortexRecord record = records.get(kmer);
        String edges = record.getEdges(color);

        return getNextOutKmers(kmer, edges);
    }

    private String getNextKmer(List<String> rawKmers, Set<String> referenceKmers) {
        if (rawKmers.size() == 1) {
            return rawKmers.get(0);
        }

        int numKmersSpanningFork = 0;
        String rawKmerSpanningFork = null;

        for (String rawKmer : rawKmers) {
            String fwKmer = SequenceUtils.alphanumericallyLowestOrientation(rawKmer);

            if (referenceKmers.contains(fwKmer)) {
                numKmersSpanningFork++;
                rawKmerSpanningFork = rawKmer;
            }
        }

        if (numKmersSpanningFork == 1) {
            return rawKmerSpanningFork;
        }

        return null;
    }

    private String getInKmerWithMatchingOrientation(int color, String kmer, Set<String> referenceKmers, String supernode) {
        String fw = SequenceUtils.alphanumericallyLowestOrientation(kmer);
        String rc = SequenceUtils.reverseComplement(fw);

        String currentKmer = null;
        String edges = records.get(fw).getEdges(color);

        if (fw.equalsIgnoreCase(supernode.substring(0, fw.length()))) {
            currentKmer = fw;
        } else if (rc.equalsIgnoreCase(supernode.substring(0, rc.length()))) {
            currentKmer = rc;
            edges = SequenceUtils.reverseComplement(edges);
        }

        if (currentKmer != null) {
            return getNextKmer(getNextInKmers(currentKmer, edges), referenceKmers);
        }

        return null;
    }

    private String getOutKmerWithMatchingOrientation(int color, String kmer, Set<String> referenceKmers, String supernode) {
        String fw = SequenceUtils.alphanumericallyLowestOrientation(kmer);
        String rc = SequenceUtils.reverseComplement(fw);

        String currentKmer = null;
        String edges = records.get(fw).getEdges(color);

        if (fw.equalsIgnoreCase(supernode.substring(supernode.length() - fw.length(), supernode.length()))) {
            currentKmer = fw;
        } else if (rc.equalsIgnoreCase(supernode.substring(supernode.length() - rc.length(), supernode.length()))) {
            currentKmer = rc;
            edges = SequenceUtils.reverseComplement(edges);
        }

        if (currentKmer != null) {
            return getNextKmer(getNextOutKmers(currentKmer, edges), referenceKmers);
        }

        return null;
    }

    public Collection<String> getReferenceGuidedSupernodes(int color, Set<String> referenceKmers) {
        Set<String> supernodes = new HashSet<String>();
        Set<String> seenKmers = new HashSet<String>();

        for (String fwRefKmer : referenceKmers) {
            if (records.containsKey(fwRefKmer) && !seenKmers.contains(fwRefKmer)) {
                String supernode = fwRefKmer;

                seenKmers.add(fwRefKmer);

                String inKmer = getNextKmer(getNextInKmers(color, fwRefKmer), referenceKmers);

                while (inKmer != null) {
                    supernode = inKmer.charAt(0) + supernode;

                    String fw = SequenceUtils.alphanumericallyLowestOrientation(inKmer);

                    if (records.containsKey(fw) && !seenKmers.contains(fw)) {
                        inKmer = getInKmerWithMatchingOrientation(color, inKmer, referenceKmers, supernode);
                    } else {
                        inKmer = null;
                    }

                    seenKmers.add(fw);
                }

                String outKmer = getNextKmer(getNextOutKmers(color, fwRefKmer), referenceKmers);

                while (outKmer != null) {
                    supernode = supernode + outKmer.charAt(outKmer.length() - 1);

                    String fw = SequenceUtils.alphanumericallyLowestOrientation(outKmer);

                    if (records.containsKey(fw) && !seenKmers.contains(fw)) {
                        outKmer = getOutKmerWithMatchingOrientation(color, outKmer, referenceKmers, supernode);
                    } else {
                        outKmer = null;
                    }

                    seenKmers.add(fw);
                }

                supernodes.add(SequenceUtils.alphanumericallyLowestOrientation(supernode));
            }
        }

        return supernodes;
    }

    public String getSupernode(int color, String startingKmer) {
        String superNode = null;

        if (records.containsKey(startingKmer)) {
            superNode = startingKmer;

            CortexRecord startingRecord = records.get(startingKmer);
            String startingEdges = startingRecord.getEdges()[color];

            // First do in kmers
            List<String> inRawKmers = new ArrayList<String>();
            for (int i = 0; i < 4; i++) {
                if (startingEdges.charAt(i) != '.') {
                    String inRawKmer = (startingEdges.charAt(i) + startingKmer.substring(0, startingKmer.length() - 1)).toUpperCase();

                    inRawKmers.add(inRawKmer);
                }
            }

            HashSet<String> seenKmers = new HashSet<String>();

            while (inRawKmers.size() == 1) {
                String inRawKmer = inRawKmers.get(0);
                inRawKmers.clear();

                superNode = inRawKmer.charAt(0) + superNode;

                String fw = SequenceUtils.alphanumericallyLowestOrientation(inRawKmer);

                if (records.containsKey(fw) && !seenKmers.contains(fw)) {
                    String rc = SequenceUtils.reverseComplement(fw);

                    String currentKmer = null;

                    String edges = records.get(fw).getEdges()[color];

                    if (fw.substring(0, fw.length()).equalsIgnoreCase(superNode.substring(0, fw.length()))) {
                        currentKmer = fw;
                    } else if (rc.substring(0, rc.length()).equalsIgnoreCase(superNode.substring(0, rc.length()))) {
                        currentKmer = rc;
                        edges = SequenceUtils.reverseComplement(edges);
                    }

                    if (currentKmer != null) {
                        for (int i = 0; i < 4; i++) {
                            if (edges.charAt(i) != '.') {
                                String newInRawKmer = (edges.charAt(i) + currentKmer.substring(0, currentKmer.length() - 1)).toUpperCase();
                                inRawKmers.add(newInRawKmer);
                            }
                        }
                    }

                    seenKmers.add(currentKmer);
                }
            }

            // Now do out kmers
            List<String> outRawKmers = new ArrayList<String>();
            for (int i = 4; i < 8; i++) {
                if (startingEdges.charAt(i) != '.') {
                    String outRawKmer = (startingKmer.substring(1, startingKmer.length()) + startingEdges.charAt(i)).toUpperCase();

                    outRawKmers.add(outRawKmer);
                }
            }

            seenKmers.clear();

            while (outRawKmers.size() == 1) {
                String outRawKmer = outRawKmers.get(0);
                outRawKmers.clear();

                superNode = superNode + outRawKmer.charAt(outRawKmer.length() - 1);

                String fw = SequenceUtils.alphanumericallyLowestOrientation(outRawKmer);

                if (records.containsKey(fw) && !seenKmers.contains(fw)) {
                    String rc = SequenceUtils.reverseComplement(fw);

                    String currentKmer = null;

                    String edges = records.get(fw).getEdges()[color];

                    if (fw.substring(0, fw.length()).equalsIgnoreCase(superNode.substring(superNode.length() - fw.length(), superNode.length()))) {
                        currentKmer = fw;
                    } else if (rc.substring(0, rc.length()).equalsIgnoreCase(superNode.substring(superNode.length() - rc.length(), superNode.length()))) {
                        currentKmer = rc;
                        edges = SequenceUtils.reverseComplement(edges);
                    }

                    if (currentKmer != null) {
                        for (int i = 4; i < 8; i++) {
                            if (edges.charAt(i) != '.') {
                                String newOutRawKmer = (currentKmer.substring(1, currentKmer.length()) + edges.charAt(i)).toUpperCase();
                                outRawKmers.add(newOutRawKmer);
                            }
                        }
                    }

                    seenKmers.add(currentKmer);
                }
            }
        }

        return SequenceUtils.alphanumericallyLowestOrientation(superNode);
    }
}
