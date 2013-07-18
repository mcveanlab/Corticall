package uk.ac.ox.well.indiana.analyses.kmerSharing;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.tools.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class JoinSupernodes extends Module {
    @Argument(fullName="supernodes", shortName="sn", doc="File of supernodes")
    public File SUPERNODES;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size to use for alignment")
    public Integer KMER_SIZE = 20;

    @Output
    public PrintStream out;

    private class ExtendedSupernodeInfo {
        public String extendedSupernode;
        public Set<String> usedSupernodes;
    }

    private List<String> loadSupernodes() {
        Set<String> supernodes = new HashSet<String>();

        TableReader tr = new TableReader(SUPERNODES);
        for (Map<String, String> te : tr) {
            supernodes.add(SequenceUtils.alphanumericallyLowestOrientation(te.get("superNode")));
        }

        List<String> supernodesList = new ArrayList<String>(supernodes);

        return supernodesList;
    }

    private Map<String, Set<Integer>> hashSupernodes(List<String> supernodes, int kmerSize) {
        Map<String, Set<Integer>> supernodeKmerHash = new HashMap<String, Set<Integer>>();

        for (int i = 0; i < supernodes.size(); i++) {
            String supernode = supernodes.get(i);
            for (int j = 0; j <= supernode.length() - kmerSize; j++) {
                String fw = SequenceUtils.alphanumericallyLowestOrientation(supernode.substring(j, j + kmerSize));

                if (!supernodeKmerHash.containsKey(fw)) {
                    supernodeKmerHash.put(fw, new HashSet<Integer>());
                }

                supernodeKmerHash.get(fw).add(i);
            }
        }

        return supernodeKmerHash;
    }

    private Map<String, String> loadUniqueKmers() {
        Map<String, String> uniqueKmers = new HashMap<String, String>();

        TableReader tr = new TableReader(SUPERNODES);
        for (Map<String, String> te : tr) {
            String kmer = te.get("kmer");
            String gene = te.get("genes");

            uniqueKmers.put(kmer, gene);
        }

        return uniqueKmers;
    }

    private int longestContigLength(Collection<String> extendedSupernodes) {
        int length = 0;

        for (String extendedSupernode : extendedSupernodes) {
            if (extendedSupernode.length() > length) {
                length = extendedSupernode.length();
            }
        }

        return length;
    }

    private String joinSupernodes(String s1, String s2, String seedFw) {
        String seedRc = SequenceUtils.reverseComplement(seedFw);
        String seed = seedFw;

        if ((s1.contains(seedFw) && s2.contains(seedRc)) || (s1.contains(seedRc) && s2.contains(seedFw))) {
            s2 = SequenceUtils.reverseComplement(s2);
        }

        if (s1.contains(seedRc)) {
            seed = seedRc;
        }

        String s3 = s1;

        int seedPos2 = s2.indexOf(seed);

        if (seedPos2 > -1) {
            StringBuilder overlapBuilder1 = new StringBuilder();
            StringBuilder overlapBuilder2 = new StringBuilder();

            for (int i = 0, j = seedPos2; i < s1.length() && j < s2.length(); i++, j++) {
                overlapBuilder1.append(s1.charAt(i));
                overlapBuilder2.append(s2.charAt(j));
            }

            String overlap1 = overlapBuilder1.toString();
            String overlap2 = overlapBuilder2.toString();


            if (overlap1.equals(overlap2)) {
                s2 = s2.replaceAll(overlap1, "");

                s3 = s2 + s1;
            }
        }

        return s3;
    }

    private String lookForSupernodeExtensions(String supernode, List<String> supernodes, Map<String, Set<Integer>> supernodeKmerHash, Set<String> seenSupernodes) {
        String firstKmerFw = SequenceUtils.alphanumericallyLowestOrientation(supernode.substring(0, KMER_SIZE));
        Set<Integer> firstKmerCandidates = supernodeKmerHash.get(firstKmerFw);

        String joinedSupernodes = supernode;

        for (Integer firstKmerCandidate : firstKmerCandidates) {
            String candidate = supernodes.get(firstKmerCandidate);

            if (!joinedSupernodes.equals(candidate) && !seenSupernodes.contains(candidate)) {
                joinedSupernodes = joinSupernodes(joinedSupernodes, candidate, firstKmerFw);

                seenSupernodes.add(candidate);
            }
        }

        return joinedSupernodes;
    }

    private String extendSupernode(String supernode, List<String> supernodes, Map<String, Set<Integer>> supernodeKmerHash, Set<String> seenSupernodes) {
        String joinedSupernodes = lookForSupernodeExtensions(supernode, supernodes, supernodeKmerHash, seenSupernodes);

        return joinedSupernodes;
    }

    @Override
    public void execute() {
        List<String> supernodes = loadSupernodes();
        Set<String> extendedSupernodes = new HashSet<String>();

        while (supernodes.size() > extendedSupernodes.size()) {
            if (extendedSupernodes.size() != 0) {
                supernodes.clear();
                supernodes.addAll(extendedSupernodes);
                extendedSupernodes.clear();
            }

            log.info("num_supernodes={} n50_length={} n50_value={}", supernodes.size(), SequenceUtils.computeN50Length(supernodes), SequenceUtils.computeN50Value(supernodes));

            Map<String, Set<Integer>> supernodeKmerHash = hashSupernodes(supernodes, KMER_SIZE);
            Set<String> seenSupernodes = new HashSet<String>();

            for (String supernode : supernodes) {
                if (!seenSupernodes.contains(supernode)) {
                    seenSupernodes.add(supernode);

                    log.debug("supernode: {}", supernode);

                    String oldSupernode = supernode;

                    boolean stillExtending = true;

                    while (stillExtending) {
                        String newSupernode = extendSupernode(oldSupernode, supernodes, supernodeKmerHash, seenSupernodes);

                        log.debug(" extended: {}", newSupernode);

                        oldSupernode = newSupernode;

                        if (newSupernode.equals(oldSupernode)) {
                            stillExtending = false;
                        }
                    }

                    oldSupernode = SequenceUtils.reverseComplement(oldSupernode);

                    stillExtending = true;

                    while (stillExtending) {
                        String newSupernode = extendSupernode(oldSupernode, supernodes, supernodeKmerHash, seenSupernodes);

                        log.debug(" extended: {}", newSupernode);

                        oldSupernode = newSupernode;

                        if (newSupernode.equals(oldSupernode)) {
                            stillExtending = false;
                        }
                    }

                    log.debug("    final: {}", oldSupernode);
                    log.debug("");

                    extendedSupernodes.add(oldSupernode);
                }
            }
        }

        Map<String, String> uniqueKmers = loadUniqueKmers();
        int kmerSize = uniqueKmers.keySet().iterator().next().length();

        out.println("color\tkmer\tgenes\tdomains\tsuperNode");

        for (String extendedSupernode : extendedSupernodes) {
            List<String> kmers = new ArrayList<String>();
            List<String> genes = new ArrayList<String>();

            for (int i = 0; i <= extendedSupernode.length() - kmerSize; i++) {
                String fw = SequenceUtils.alphanumericallyLowestOrientation(extendedSupernode.substring(i, i + kmerSize));

                if (uniqueKmers.containsKey(fw)) {
                    kmers.add(fw);
                    genes.add(uniqueKmers.get(fw));
                }
            }

            String kmersString = Joiner.on(",").join(kmers);
            String genesString = Joiner.on(",").join(genes);

            out.printf("%d\t%s\t%s\tunknown\t%s\n", 0, kmersString, genesString, extendedSupernode);
        }
    }
}
