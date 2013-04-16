package uk.ac.ox.well.indiana.analyses.kmerSharing;

import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.utils.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class FindMeganode extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="stoppingThreshold", shortName="st", doc="Number of variants to allow before we stop")
    public Integer STOPPING_THRESHOLD = 1;

    @Argument(fullName="kmerReferencePanel", shortName="krp", doc="Kmer reference panel")
    public File KMER_REFERENCE_PANEL;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        HashMap<String, CortexRecord> kmers = new HashMap<String, CortexRecord>();

        int numRecords = 0;
        for (CortexRecord cr : CORTEX_GRAPH) {
            if (numRecords % (CORTEX_GRAPH.getNumRecords() / 10) == 0) {
                log.info("Processed {}/{} records", numRecords, CORTEX_GRAPH.getNumRecords());
            }
            numRecords++;

            kmers.put(cr.getKmerString(), cr);
        }

        TableReader table = new TableReader(KMER_REFERENCE_PANEL);
        for (HashMap<String, String> entry : table) {
            String kmer = entry.get("kmer");

            log.info("kmer: {}", kmer);

            String supernode = getSupernode(kmer, 0, kmers);
            log.info("supernode: {}", supernode);

            getMeganode(0, kmer, kmers, STOPPING_THRESHOLD);
            log.info("");
        }
    }

    public void getMeganode(int color, String startingKmer, HashMap<String, CortexRecord> records, int stoppingThreshold) {
        HashSet<String> seenKmers = new HashSet<String>();

        Set<String> branches = walkGraphToLeft(color, startingKmer, records, stoppingThreshold, seenKmers);

        for (String branch : branches) {
            log.info("branch: {}", branch);
        }
    }

    private Set<String> walkGraphToLeft(int color, String supernode, HashMap<String, CortexRecord> records, int stoppingThreshold, HashSet<String> seenKmers) {
        Set<String> branches = new HashSet<String>();

        String kmer = supernode.substring(0, CORTEX_GRAPH.getKmerSize());
        String fw = SequenceUtils.getAlphanumericallyLowestOrientation(kmer);

        while (fw != null && !seenKmers.contains(fw) && records.containsKey(fw)) {
            String rc = SequenceUtils.getReverseComplement(fw);

            CortexRecord cr = records.get(fw);
            String edges = cr.getEdges()[color];

            if (rc.equalsIgnoreCase(kmer)) {
                edges = SequenceUtils.getReverseComplement(edges);
            }

            List<String> inKmers = new ArrayList<String>();

            for (int i = 0; i < 4; i++) {
                if (edges.charAt(i) != '.') {
                    String newKmer = (edges.charAt(i) + kmer.substring(0, kmer.length() - 1)).toUpperCase();

                    inKmers.add(newKmer);
                }
            }

            seenKmers.add(fw);

            if (inKmers.size() == 1) {
                String inKmer = inKmers.get(0);
                supernode = inKmer.charAt(0) + supernode;

                kmer = supernode.substring(0, CORTEX_GRAPH.getKmerSize());
                fw = SequenceUtils.getAlphanumericallyLowestOrientation(kmer);
            } else {
                fw = null;

                //branches.add(supernode);

                if (stoppingThreshold > 0) {
                    for (String inKmer : inKmers) {
                        //Set<String> subBranches = walkGraphToLeft(color, inKmer, records, stoppingThreshold - 1, seenKmers);
                        Set<String> subBranches = walkGraphToLeft(color, inKmer, records, stoppingThreshold - 1, new HashSet<String>());
                        branches.addAll(subBranches);
                    }
                }
            }
        }

        branches.add(supernode);

        return branches;
    }

    public void walkGraphToRight(String supernode, HashMap<String, CortexRecord> records, int stoppingThreshold) {

    }

    public String getSupernode(String startingKmer, int color, HashMap<String, CortexRecord> records) {
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

                String fw = SequenceUtils.getAlphanumericallyLowestOrientation(inRawKmer);

                if (records.containsKey(fw) && !seenKmers.contains(fw)) {
                    String rc = SequenceUtils.getReverseComplement(fw);

                    String currentKmer = null;

                    String edges = records.get(fw).getEdges()[color];

                    if (fw.substring(0, fw.length()).equalsIgnoreCase(superNode.substring(0, fw.length()))) {
                        currentKmer = fw;
                    } else if (rc.substring(0, rc.length()).equalsIgnoreCase(superNode.substring(0, rc.length()))) {
                        currentKmer = rc;
                        edges = SequenceUtils.getReverseComplement(edges);
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

                String fw = SequenceUtils.getAlphanumericallyLowestOrientation(outRawKmer);

                if (records.containsKey(fw) && !seenKmers.contains(fw)) {
                    String rc = SequenceUtils.getReverseComplement(fw);

                    String currentKmer = null;

                    String edges = records.get(fw).getEdges()[color];

                    if (fw.substring(0, fw.length()).equalsIgnoreCase(superNode.substring(superNode.length() - fw.length(), superNode.length()))) {
                        currentKmer = fw;
                    } else if (rc.substring(0, rc.length()).equalsIgnoreCase(superNode.substring(superNode.length() - rc.length(), superNode.length()))) {
                        currentKmer = rc;
                        edges = SequenceUtils.getReverseComplement(edges);
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

        return superNode;
    }
}
