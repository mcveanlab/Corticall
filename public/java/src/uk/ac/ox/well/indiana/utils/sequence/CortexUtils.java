package uk.ac.ox.well.indiana.utils.sequence;

import htsjdk.samtools.util.SequenceUtil;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksMap;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;

import java.util.*;

/**
 * A set of utilities for dealing with Cortex graphs and records
 */
public class CortexUtils {
    /**
     * Private constructor - this class cannot be instantiated!
     */
    private CortexUtils() {}

    /**
     * Return the next kmer, in the orientation of the given kmer
     *
     * @param cg    the sorted Cortex graph
     * @param kmer  the kmer
     * @return      the next kmer (in the orientation of the given kmer)
     */
    public static String getNextKmer(CortexGraph cg, String kmer) {
        CortexKmer ck = new CortexKmer(kmer);
        CortexRecord cr = cg.findRecord(ck);

        if (cr != null) {
            Collection<String> outEdges = cr.getOutEdgesAsStrings(0);

            if (ck.isFlipped()) {
                outEdges = cr.getInEdgesComplementAsStrings(0);
            }

            if (outEdges.size() == 1) {
                String outEdge = outEdges.iterator().next();

                return kmer.substring(1, kmer.length()) + outEdge;
            }
        }

        return null;
    }

    /**
     * Return the next kmers, in the orientation of the given kmer
     *
     * @param cg    the sorted Cortex graph
     * @param kmer  the kmer
     * @return      the next kmers (in the orientation of the given kmer)
     */
    public static Set<String> getNextKmers(CortexGraph cg, String kmer) {
        CortexKmer ck = new CortexKmer(kmer);
        CortexRecord cr = cg.findRecord(ck);
        Set<String> nextKmers = new HashSet<String>();

        if (cr != null) {
            Collection<String> outEdges = cr.getOutEdgesAsStrings(0);

            if (ck.isFlipped()) {
                outEdges = cr.getInEdgesComplementAsStrings(0);
            }


            for (String outEdge : outEdges) {
                nextKmers.add(kmer.substring(1, kmer.length()) + outEdge);
            }
        }

        return nextKmers;
    }

    /**
     * Return the previous kmer, in the orientation of the given kmer
     *
     * @param cg    the sorted Cortex graph
     * @param kmer  the kmer
     * @return      the preivous kmer (in the orientation of the given kmer)
     */
    public static String getPrevKmer(CortexGraph cg, String kmer) {
        CortexKmer ck = new CortexKmer(kmer);
        CortexRecord cr = cg.findRecord(ck);

        if (cr != null) {
            Collection<String> inEdges = cr.getInEdgesAsStrings(0);

            if (ck.isFlipped()) {
                inEdges = cr.getOutEdgesComplementAsStrings(0);
            }

            if (inEdges.size() == 1) {
                String inEdge = inEdges.iterator().next();

                return inEdge + kmer.substring(0, kmer.length() - 1);
            }
        }

        return null;
    }

    /**
     * Return the previous kmers, in the orientation of the given kmer
     *
     * @param cg    the sorted Cortex graph
     * @param kmer  the kmer
     * @return      the preivous kmers (in the orientation of the given kmer)
     */
    public static Set<String> getPrevKmers(CortexGraph cg, String kmer) {
        CortexKmer ck = new CortexKmer(kmer);
        CortexRecord cr = cg.findRecord(ck);
        Set<String> prevKmers = new HashSet<String>();

        if (cr != null) {
            Collection<String> inEdges = cr.getInEdgesAsStrings(0);

            if (ck.isFlipped()) {
                inEdges = cr.getOutEdgesComplementAsStrings(0);
            }

            for (String inEdge : inEdges) {
                prevKmers.add(inEdge + kmer.substring(0, kmer.length() - 1));
            }
        }

        return prevKmers;
    }

    /**
     * Flip the information in a links record to the opposite orientation
     *
     * @param clr  the Cortex links record
     * @return     the flipped Cortex links record
     */
    public static CortexLinksRecord flipLinksRecord(CortexLinksRecord clr) {
        String rsk = SequenceUtils.reverseComplement(clr.getKmerAsString());
        List<CortexJunctionsRecord> cjrs = new ArrayList<CortexJunctionsRecord>();
        for (CortexJunctionsRecord cjr : clr.getJunctions()) {
            cjrs.add(flipJunctionsRecord(cjr));
        }

        return new CortexLinksRecord(rsk, cjrs);
    }

    /**
     * Flip the information in a junctions record to the opposite orientation
     *
     * @param cjr  the Cortex link junctions record
     * @return     the flipped Cortex link junctions record
     */
    public static CortexJunctionsRecord flipJunctionsRecord(CortexJunctionsRecord cjr) {
        return new CortexJunctionsRecord(!cjr.isForward(),
                                         cjr.getNumKmers(),
                                         cjr.getNumJunctions(),
                                         cjr.getCoverages(),
                                         SequenceUtils.complement(cjr.getJunctions()),
                                         SequenceUtils.reverseComplement(cjr.getSeq()),
                                         cjr.getJunctionPositions()
                   );
    }

    /**
     * Adjust the information in a links record to the opposite orientation
     *
     * @param sk   the kmer in desired orientation
     * @param clr  the Cortex link junctions record
     * @return     the reoriented Cortex links record
     */
    public static CortexLinksRecord orientLinksRecord(String sk, CortexLinksRecord clr) {
        List<CortexJunctionsRecord> cjrs = new ArrayList<CortexJunctionsRecord>();
        for (CortexJunctionsRecord cjr : clr.getJunctions()) {
            cjrs.add(orientJunctionsRecord(sk, cjr));
        }

        return new CortexLinksRecord(sk, cjrs);
    }

    /**
     * Adjust the information in a junctions to facilitate navigation in contig orientation
     *
     * @param sk   the kmer in desired orientation
     * @param cjr  the Cortex link junctions record
     * @return     the reoriented Cortex link junctions record
     */
    public static CortexJunctionsRecord orientJunctionsRecord(String sk, CortexJunctionsRecord cjr) {
        CortexKmer ck = new CortexKmer(sk);
        boolean isForward = cjr.isForward();
        String junctions = cjr.getJunctions();

        if (!ck.isFlipped() && isForward) { // --> F
            // do nothing
        } else if (!ck.isFlipped() && !isForward) { // --> R
            junctions = SequenceUtils.complement(junctions);
        } else if (ck.isFlipped() && isForward) { // <-- F
            isForward = !isForward;
            junctions = SequenceUtils.complement(junctions);
        } else if (ck.isFlipped() && !isForward) { // <-- R
            isForward = !isForward;
        }

        return new CortexJunctionsRecord(isForward, cjr.getNumKmers(), cjr.getNumJunctions(), cjr.getCoverages(), junctions, cjr.getSeq(), cjr.getJunctionPositions());
    }

    /**
     * Return a list of all of the kmers in the graph between a starting kmer and the end of the junctions record
     * @param cg   the Cortex graph
     * @param sk   the kmer in desired orientation
     * @param cjr  the Cortex link junctions record
     * @return     the list of kmers
     */
    public static List<String> getKmersInLinkByNavigation(CortexGraph cg, String sk, CortexJunctionsRecord cjr) {
        cjr = orientJunctionsRecord(sk, cjr);

        int junctionsUsed = 0;
        String curKmer = sk;
        boolean isForward = cjr.isForward();
        String junctions = cjr.getJunctions();

        List<String> kmersInLink = new ArrayList<String>();
        kmersInLink.add(sk);

        if (isForward) {
            Set<String> nextKmers = CortexUtils.getNextKmers(cg, sk);

            while (nextKmers.size() > 0 && junctionsUsed < junctions.length()) {
                if (nextKmers.size() == 1) {
                    curKmer = nextKmers.iterator().next();
                    nextKmers = CortexUtils.getNextKmers(cg, curKmer);

                    kmersInLink.add(curKmer);
                } else {
                    char nbase = junctions.charAt(junctionsUsed);

                    String expectedNextKmer = curKmer.substring(1, curKmer.length()) + nbase;

                    if (nextKmers.contains(expectedNextKmer)) {
                        curKmer = expectedNextKmer;
                        nextKmers = CortexUtils.getNextKmers(cg, curKmer);

                        kmersInLink.add(expectedNextKmer);
                    } else {
                        throw new IndianaException("Junction record specified a navigation that conflicted with the graph: " + junctions + " " + nbase + " " + expectedNextKmer + " " + nextKmers);
                    }

                    junctionsUsed++;
                }
            }
        } else {
            Set<String> prevKmers = CortexUtils.getPrevKmers(cg, sk);

            while (prevKmers.size() > 0 && junctionsUsed < junctions.length()) {
                if (prevKmers.size() == 1) {
                    curKmer = prevKmers.iterator().next();
                    prevKmers = CortexUtils.getPrevKmers(cg, curKmer);

                    kmersInLink.add(0, curKmer);
                } else {
                    char pbase = junctions.charAt(junctionsUsed);

                    String expectedPrevKmer = pbase + curKmer.substring(0, curKmer.length() - 1);

                    if (prevKmers.contains(expectedPrevKmer)) {
                        curKmer = expectedPrevKmer;
                        prevKmers = CortexUtils.getPrevKmers(cg, curKmer);

                        kmersInLink.add(0, expectedPrevKmer);
                    } else {
                        throw new IndianaException("Junction record specified a navigation that conflicted with the graph: " + junctions + " " + pbase + " " + expectedPrevKmer + " " + prevKmers);
                    }

                    junctionsUsed++;
                }
            }
        }

        return kmersInLink;
    }

    public static List<String> getKmersInLinkFromSeq(CortexGraph cg, String sk, CortexJunctionsRecord cjr) {
        List<String> expectedKmers = new ArrayList<String>();
        String seq = cjr.getSeq();
        CortexKmer ck = new CortexKmer(sk);
        if ((!ck.isFlipped() && !cjr.isForward()) || (ck.isFlipped() && cjr.isForward())) {
            seq = SequenceUtils.reverseComplement(seq);
        }
        for (int i = 0; i <= seq.length() - cg.getKmerSize(); i++) {
            String kmer = seq.substring(i, i + cg.getKmerSize());
            expectedKmers.add(kmer);
        }

        return expectedKmers;
    }

    public static List<String> getKmersInLink(CortexGraph cg, String sk, CortexJunctionsRecord cjr) {
        return (cjr.getSeq() != null) ? getKmersInLinkFromSeq(cg, sk, cjr) : getKmersInLinkByNavigation(cg, sk, cjr);
    }

    public static Map<String, Integer> getKmersAndCoverageInLink(CortexGraph cg, String sk, CortexLinksRecord clr) {
        Map<String, Integer> kmersAndCoverage = new HashMap<String, Integer>();

        for (CortexJunctionsRecord cjr : clr.getJunctions()) {
            List<String> kil = getKmersInLink(cg, sk, cjr);

            for (String ki : kil) {
                if (!kmersAndCoverage.containsKey(ki)) {
                    kmersAndCoverage.put(ki, cjr.getCoverage(0));
                } else {
                    kmersAndCoverage.put(ki, cjr.getCoverage(0) + kmersAndCoverage.get(ki));
                }
            }
        }

        return kmersAndCoverage;
    }
}
