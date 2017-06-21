package uk.ac.ox.well.indiana.commands.caller.call;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.jgrapht.graph.EdgeReversedGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksMap;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.stoppingconditions.ContigStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.NahrStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.TraversalStopper;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

/**
 * Created by kiran on 20/06/2017.
 */
public class CallNahrEvents extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Links")
    public CortexLinksMap LINKS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="refs", shortName="R", doc="References")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(PARENTS);
        List<Integer> recruitColors = GRAPH.getColorsForSampleNames(new ArrayList<>(LOOKUPS.keySet()));

        logColorAssignments(childColor, parentColors, recruitColors);

        Map<CortexKmer, Boolean> usedRois = loadRois();

        TraversalEngine ce = initializeTraversalEngine(childColor, parentColors, recruitColors, ContigStopper.class);
        //TraversalEngine ne = initializeTraversalEngine(childColor, parentColors, recruitColors, NahrStopper.class);

        Map<String, String> contigEncoding = new HashMap<>();
        String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
        Random r = new Random();
        for (String key : LOOKUPS.keySet()) {
            Set<String> usedCodes = new HashSet<>();

            for (SAMSequenceRecord ssr : LOOKUPS.get(key).getReferenceSequence().getSequenceDictionary().getSequences()) {
                String sn = ssr.getSequenceName();
                String c;
                do {
                    c = String.valueOf(alphabet.charAt(r.nextInt(alphabet.length())));
                } while(usedCodes.contains(c) && usedCodes.size() < alphabet.length());

                if (!usedCodes.contains(c)) {
                    contigEncoding.put(sn, c);
                    usedCodes.add(c);
                }
            }
        }

        //Pattern motif = Pattern.compile("\\.+_*(\\w)\\1+");
        Pattern motif = Pattern.compile("\\.+_*(\\w)\1+");

        for (CortexRecord rr : ROI) {
            if (!usedRois.get(rr.getCortexKmer())) {
                DirectedWeightedPseudograph<CortexVertex, CortexEdge> cg = ce.dfs(rr.getKmerAsString());
                CortexVertex rv = new CortexVertex(rr.getKmerAsString(), GRAPH.findRecord(rr.getCortexKmer()));

                for (String key : LOOKUPS.keySet()) {
                    KmerLookup kl = LOOKUPS.get(key);

                    DepthFirstIterator<CortexVertex, CortexEdge> dfsf = new DepthFirstIterator<>(cg, rv);
                    String fContigCount = getContigCounts(kl, dfsf, usedRois, contigEncoding, true);

                    DepthFirstIterator<CortexVertex, CortexEdge> dfsr = new DepthFirstIterator<>(new EdgeReversedGraph<>(cg), rv);
                    String rContigCount = getContigCounts(kl, dfsr, usedRois, contigEncoding, true);

                    log.info("{} rContigCount: {} {}", rr.getCortexKmer(), key, rContigCount);
                    log.info("{} fContigCount: {} {}", rr.getCortexKmer(), key, fContigCount);

                    Matcher rMatcher = motif.matcher(rContigCount);
                    log.info("{} {}", rMatcher, rMatcher.groupCount());

                    Matcher fMatcher = motif.matcher(fContigCount);
                    log.info("{} {}", fMatcher, fMatcher.groupCount());

                    rContigCount = SequenceUtils.reverse(rContigCount);
                    String contig = rContigCount + fContigCount;
                    log.info("{}       contig: {} {}", rr.getCortexKmer(), key, contig);

                    log.info("");
                }

                for (CortexVertex cv : cg.vertexSet()) {
                    usedRois.put(cv.getCk(), true);
                }
            }
        }
    }

    private String getContigCounts(KmerLookup kl, DepthFirstIterator<CortexVertex, CortexEdge> dfs, Map<CortexKmer, Boolean> usedRois, Map<String, String> contigEncoding, boolean expand) {
        StringBuilder sb = new StringBuilder();

        while (dfs.hasNext()) {
            CortexVertex cv = dfs.next();
            Set<Interval> loci = kl.findKmer(cv.getSk());
            if (usedRois.containsKey(cv.getCk())) {
                if (expand) {
                    sb.append(".");
                }
            } else if (loci.size() == 1) {
                for (Interval locus : loci) {
                    sb.append(contigEncoding.get(locus.getContig()));
                }
            } else {
                if (expand) {
                    sb.append("_");
                }
            }
        }

        return sb.toString();
    }

    private TraversalEngine initializeTraversalEngine(int childColor, List<Integer> parentColors, List<Integer> recruitColors, Class<? extends TraversalStopper<CortexVertex, CortexEdge>> stoppingRule) {
        TraversalEngine e = new TraversalEngineFactory()
            .combinationOperator(AND)
            .traversalDirection(BOTH)
            .traversalColor(childColor)
            .joiningColors(parentColors)
            .recruitmentColors(recruitColors)
            .rois(ROI)
            .connectAllNeighbors(true)
            .stopper(stoppingRule)
            .graph(GRAPH)
            .make();

        return e;
    }

    private Map<CortexKmer, Boolean> loadRois() {
        Map<CortexKmer, Boolean> usedRois = new HashMap<>();
        for (CortexRecord rr : ROI) {
            usedRois.put(rr.getCortexKmer(), false);
        }

        return usedRois;
    }

    private void logColorAssignments(int childColor, List<Integer> parentColors, List<Integer> recruitColors) {
        log.info("Colors:");
        log.info("    child: {} {}", childColor, CHILD);
        for (int i = 0; i < parentColors.size(); i++) {
            log.info("   parent: {} {}", parentColors.get(i), GRAPH.getSampleName(parentColors.get(i)));
        }
        for (int i = 0; i < recruitColors.size(); i++) {
            log.info("  recruit: {} {}", recruitColors.get(i), GRAPH.getSampleName(recruitColors.get(i)));
        }
    }
}
