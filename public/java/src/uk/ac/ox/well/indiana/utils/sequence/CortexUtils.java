package uk.ac.ox.well.indiana.utils.sequence;

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
                                         SequenceUtils.complement(cjr.getJunctions()));
    }

    /**
     * Return a list of all of the kmers in the graph between a starting kmer and the end of the junctions record
     * @param cg   the Cortex graph
     * @param sk   the kmer in desired orientation
     * @param cjr  the Cortex link junctions record
     * @return     the list of kmers
     */
    public static List<String> getKmersInLink(CortexGraph cg, String sk, CortexJunctionsRecord cjr) {
        CortexKmer ck = new CortexKmer(sk);

        if (ck.isFlipped()) {
            cjr = CortexUtils.flipJunctionsRecord(cjr);
        }

        String junctions = cjr.getJunctions();
        int junctionsUsed = 0;
        String curKmer = sk;

        List<String> kmersInLink = new ArrayList<String>();
        kmersInLink.add(sk);

        if (cjr.isForward()) {
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
                        throw new IndianaException("Junction record specified a navigation that conflicted with the graph: " + junctions + " " + nbase);
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

                    //kmersInLink.add(curKmer);
                    kmersInLink.add(0, curKmer);
                } else {
                    char pbase = junctions.charAt(junctionsUsed);

                    String expectedPrevKmer = pbase + curKmer.substring(0, curKmer.length() - 1);

                    if (prevKmers.contains(expectedPrevKmer)) {
                        curKmer = expectedPrevKmer;
                        prevKmers = CortexUtils.getPrevKmers(cg, curKmer);

                        //kmersInLink.add(expectedPrevKmer);
                        kmersInLink.add(0, expectedPrevKmer);
                    } else {
                        throw new IndianaException("Junction record specified a navigation that conflicted with the graph: " + junctions + " " + pbase);
                    }

                    junctionsUsed++;
                }
            }
        }

        return kmersInLink;
    }

}
