package uk.ac.ox.well.cortexjdk.utils.traversal;

import org.jgrapht.Graphs;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.util.*;

public class TraversalUtils {
    private TraversalUtils() {}

    public static DirectedWeightedPseudograph<CortexVertex, CortexEdge> toGraph(List<CortexVertex> walk, int color) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = new DirectedWeightedPseudograph<>(CortexEdge.class);

        CortexVertex pv = walk.get(0);
        g.addVertex(pv);

        for (int i = 1; i < walk.size(); i++) {
            CortexVertex nv = walk.get(i);

            g.addVertex(nv);
            g.addEdge(pv, nv, new CortexEdge(pv, nv, color, 1.0));
        }

        return g;
    }

    public static String toContig(List<CortexVertex> walk) {
        StringBuilder sb = new StringBuilder();

        for (CortexVertex cv : walk) {
            String sk = cv.getKmerAsString();

            if (sb.length() == 0) {
                sb.append(sk);
            } else {
                sb.append(sk.substring(sk.length() - 1, sk.length()));
            }
        }

        return sb.toString();
    }

    public static List<CortexVertex> toWalk(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, CanonicalKmer ck) {
        return TraversalUtils.toWalk(g, ck.getKmerAsString());
    }

    public static List<CortexVertex> toWalk(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, String sk) {
        List<CortexVertex> w = new ArrayList<>();

        if (g == null) { return w; }

        CortexVertex seed = null;
        for (CortexVertex v : g.vertexSet()) {
            if (v.getKmerAsString().equals(sk) && v.getCopyIndex() == 0) {
                seed = v;
                break;
            }
        }

        if (seed != null) {
            w.add(seed);

            Set<CortexVertex> seen = new HashSet<>();
            CortexVertex cv = seed;
            while (cv != null && !seen.contains(cv)) {
                List<CortexVertex> nvs = Graphs.successorListOf(g, cv);
                nvs.remove(cv);

                CortexVertex nv = null;

                if (nvs.size() == 1) {
                    nv = nvs.get(0);
                } else if (nvs.size() > 1) {
                    boolean allKmersTheSame = true;
                    for (int i = 1; i < nvs.size(); i++) {
                        if (!nvs.get(0).getCanonicalKmer().equals(nvs.get(i).getCanonicalKmer())) {
                            allKmersTheSame = false;
                            break;
                        }
                    }

                    if (allKmersTheSame) {
                        nvs.sort((o1, o2) -> o1.getCopyIndex() < o2.getCopyIndex() ? -1 : 1);

                        if (nvs.size() > 0) {
                            nv = nvs.get(0);
                        }
                    }
                }

                if (nv != null) {
                    w.add(nv);
                    seen.add(cv);
                }

                cv = nv;
            }

            seen = new HashSet<>();
            cv = seed;
            while (cv != null && !seen.contains(cv)) {
                List<CortexVertex> pvs = Graphs.predecessorListOf(g, cv);
                pvs.remove(cv);

                CortexVertex pv = null;

                if (pvs.size() == 1) {
                    pv = pvs.get(0);
                } else if (pvs.size() > 1) {
                    boolean allKmersTheSame = true;
                    for (int i = 1; i < pvs.size(); i++) {
                        if (!pvs.get(0).getCanonicalKmer().equals(pvs.get(i).getCanonicalKmer())) {
                            allKmersTheSame = false;
                            break;
                        }
                    }

                    if (allKmersTheSame) {
                        pvs.sort((o1, o2) -> o1.getCopyIndex() > o2.getCopyIndex() ? -1 : 1);

                        if (pvs.size() > 0) {
                            pv = pvs.get(0);
                        }
                    }
                }

                if (pv != null) {
                    w.add(0, pv);
                    seen.add(cv);
                }

                cv = pv;
            }
        }

        return w;
    }

    public static CortexVertex findVertex(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, CanonicalKmer ck) {
        for (CortexVertex v : g.vertexSet()) {
            if (v.getCanonicalKmer().equals(ck)) {
                return v;
            }
        }

        return null;
    }

    public static CortexVertex findVertex(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, String sk) {
        for (CortexVertex v : g.vertexSet()) {
            if (v.getKmerAsString().equals(sk)) {
                return v;
            }
        }

        return null;
    }

    public static Map<Integer, Set<CortexByteKmer>> getAllPrevKmers(CortexRecord cr, boolean isFlipped) {
        Map<Integer, Set<CortexByteKmer>> prevKmers = new HashMap<>();

        if (cr != null) {
            byte[] sk = !isFlipped ? cr.getKmerAsBytes() : SequenceUtils.reverseComplement(cr.getKmerAsBytes());
            Map<Integer, Set<Byte>> inEdges = getInEdges(cr, isFlipped);

            for (int c = 0; c < cr.getNumColors(); c++) {
                Set<CortexByteKmer> inKmers = new HashSet<>();

                for (byte inEdge : inEdges.get(c)) {
                    byte[] inKmer = new byte[sk.length];
                    inKmer[0] = inEdge;
                    System.arraycopy(sk, 0, inKmer, 1, sk.length - 1);

                    inKmers.add(new CortexByteKmer(inKmer));
                }

                prevKmers.put(c, inKmers);
            }
        }

        return prevKmers;
    }

    public static Map<Integer, Set<CortexByteKmer>> getAllNextKmers(CortexRecord cr, boolean isFlipped) {
        Map<Integer, Set<CortexByteKmer>> nextKmers = new HashMap<>();

        if (cr != null) {
            byte[] sk = !isFlipped ? cr.getKmerAsBytes() : SequenceUtils.reverseComplement(cr.getKmerAsBytes());
            Map<Integer, Set<Byte>> outEdges = getOutEdges(cr, isFlipped);

            for (int c = 0; c < cr.getNumColors(); c++) {
                Set<CortexByteKmer> outKmers = new HashSet<>();

                for (byte outEdge : outEdges.get(c)) {
                    byte[] outKmer = new byte[sk.length];
                    System.arraycopy(sk, 1, outKmer, 0, sk.length - 1);
                    outKmer[outKmer.length - 1] = outEdge;

                    outKmers.add(new CortexByteKmer(outKmer));
                }

                nextKmers.put(c, outKmers);
            }
        }

        return nextKmers;
    }

    public static Map<Integer, Set<Byte>> getInEdges(CortexRecord cr, boolean kmerIsFlipped) {
        Map<Integer, Set<Byte>> inEdges = new HashMap<>();

        if (!kmerIsFlipped) {
            for (int c = 0; c < cr.getNumColors(); c++) {
                inEdges.put(c, new HashSet<>(cr.getInEdgesAsBytes(c, false)));
            }
        } else {
            for (int c = 0; c < cr.getNumColors(); c++) {
                inEdges.put(c, new HashSet<>(cr.getOutEdgesAsBytes(c, true)));
            }
        }

        return inEdges;
    }

    public static Map<Integer, Set<Byte>> getOutEdges(CortexRecord cr, boolean kmerIsFlipped) {
        Map<Integer, Set<Byte>> outEdges = new HashMap<>();

        if (!kmerIsFlipped) {
            for (int c = 0; c < cr.getNumColors(); c++) {
                outEdges.put(c, new HashSet<>(cr.getOutEdgesAsBytes(c, false)));
            }
        } else {
            for (int c = 0; c < cr.getNumColors(); c++) {
                outEdges.put(c, new HashSet<>(cr.getInEdgesAsBytes(c, true)));
            }
        }

        return outEdges;
    }

    public static int outDegree(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, CortexVertex v) {
        Set<CortexVertex> vs = new HashSet<>();

        Set<CortexEdge> es = g.outgoingEdgesOf(v);
        for (CortexEdge e : es) {
            vs.add(g.getEdgeTarget(e));
        }

        return vs.size();
    }

    public static int inDegree(DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, CortexVertex v) {
        Set<CortexVertex> vs = new HashSet<>();

        Set<CortexEdge> es = g.incomingEdgesOf(v);
        for (CortexEdge e : es) {
            vs.add(g.getEdgeSource(e));
        }

        return vs.size();
    }
}
