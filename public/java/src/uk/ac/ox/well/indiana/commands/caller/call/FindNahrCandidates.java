package uk.ac.ox.well.indiana.commands.caller.call;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.jetbrains.annotations.NotNull;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.jgrapht.graph.EdgeReversedGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.stoppingconditions.ContigStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.TraversalStopper;
import uk.ac.ox.well.indiana.utils.traversal.CortexEdge;
import uk.ac.ox.well.indiana.utils.traversal.CortexVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;

import java.io.PrintStream;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.AND;
import static uk.ac.ox.well.indiana.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

/**
 * Created by kiran on 20/06/2017.
 */
public class FindNahrCandidates extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="refs", shortName="R", doc="References")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Argument(fullName="nahr", shortName="n", doc="NAHR")
    public FastaSequenceFile NAHR;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(PARENTS);
        List<Integer> recruitColors = GRAPH.getColorsForSampleNames(new ArrayList<>(LOOKUPS.keySet()));

        logColorAssignments(childColor, parentColors, recruitColors);

        Map<CortexKmer, Boolean> roisThatShouldBeFound = new HashMap<>();

        ReferenceSequence rseq;
        while ((rseq = NAHR.nextSequence()) != null) {
            String seq = rseq.getBaseString();

            for (int i = 0; i <= seq.length() - ROI.getKmerSize(); i++) {
                String sk = seq.substring(i, i + ROI.getKmerSize());
                CortexKmer ck = new CortexKmer(sk);
                CortexRecord rr = ROI.findRecord(ck);

                if (rr != null) {
                    roisThatShouldBeFound.put(ck, false);
                }
            }
        }

        Map<CortexKmer, String> candidates = findNahrCandidates(childColor, parentColors, recruitColors, roisThatShouldBeFound);

        for (CortexKmer ck : candidates.keySet()) {
            if (roisThatShouldBeFound.containsKey(ck)) {
                roisThatShouldBeFound.put(ck, true);
            }
        }

        for (CortexKmer ck : roisThatShouldBeFound.keySet()) {
            log.info("novelty ck={} {}", ck, roisThatShouldBeFound.get(ck));
        }

        int i = 0;
        for (CortexKmer ck : candidates.keySet()) {
            out.println(">kmer_" + i + "_" + candidates.get(ck));
            out.println(ck.getKmerAsString());

            i++;
        }
    }

    private Map<CortexKmer, String> findNahrCandidates(int childColor, List<Integer> parentColors, List<Integer> recruitColors, Map<CortexKmer, Boolean> roisThatShouldBeFound) {
        Map<CortexKmer, String> candidates = new HashMap<>();

        String pattern = "^(\\.+)_*(([A-Za-z0-9])\\3+).*";
        Pattern motif = Pattern.compile(pattern);

        Map<String, String> contigEncoding = createContigEncoding();
        Map<CortexKmer, Boolean> usedRois = loadRois();

        TraversalEngine ce = initializeTraversalEngine(childColor, parentColors, recruitColors, ContigStopper.class);

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Finding NAHR candidates")
                .message("novel kmers processed")
                .maxRecord(ROI.getNumRecords())
                .make(log);

        for (CortexRecord rr : ROI) {
            if (!usedRois.get(rr.getCortexKmer())) {
                DirectedWeightedPseudograph<CortexVertex, CortexEdge> cg = ce.dfs(rr.getKmerAsString());
                CortexVertex rv = new CortexVertex(rr.getKmerAsString(), GRAPH.findRecord(rr.getCortexKmer()));

                for (String key : LOOKUPS.keySet()) {
                    KmerLookup kl = LOOKUPS.get(key);

                    DepthFirstIterator<CortexVertex, CortexEdge> dfsf = new DepthFirstIterator<>(cg, rv);
                    String fContigCount = getContigCounts(kl, dfsf, usedRois, contigEncoding);
                    Matcher fMatcher = motif.matcher(fContigCount);

                    DepthFirstIterator<CortexVertex, CortexEdge> dfsr = new DepthFirstIterator<>(new EdgeReversedGraph<>(cg), rv);
                    String rContigCount = getContigCounts(kl, dfsr, usedRois, contigEncoding);
                    Matcher rMatcher = motif.matcher(rContigCount);

                    boolean hasOneOfThoseKmers = false;
                    for (CortexVertex cv : cg.vertexSet()) {
                        if (roisThatShouldBeFound.containsKey(cv.getCk())) {
                            hasOneOfThoseKmers = true;
                        }
                    }

                    if (rMatcher.matches() && fMatcher.matches() &&
                        !rMatcher.group(3).equals(fMatcher.group(3)) &&
                        rMatcher.group(2).length() >= 5 && fMatcher.group(2).length() >= 5 &&
                        rMatcher.group(1).length() >= 2 && fMatcher.group(1).length() >= 2) {
                        log.info("    candidate: {} {}", rr.getCortexKmer(), key);
                        log.info("    - rContigCount: .={} {}={} {}", rMatcher.group(1).length(), rMatcher.group(3), rMatcher.group(2).length(), rContigCount);
                        log.info("    - fContigCount: .={} {}={} {}", fMatcher.group(1).length(), fMatcher.group(3), fMatcher.group(2).length(), fContigCount);

                        //candidates.put(rr.getCortexKmer(), key);
                        for (CortexVertex cv : cg.vertexSet()) {
                            if (usedRois.containsKey(cv.getCk())) {
                                candidates.put(cv.getCk(), key);
                            }
                        }
                    } else if (hasOneOfThoseKmers) {
                        log.info("    candidate: {} {}", rr.getCortexKmer(), key);
                        log.info("    - rContigCount: .={} {}={} {}", rMatcher.group(1).length(), rMatcher.group(3), rMatcher.group(2).length(), rContigCount);
                        log.info("    - fContigCount: .={} {}={} {}", fMatcher.group(1).length(), fMatcher.group(3), fMatcher.group(2).length(), fContigCount);
                    }
                }

                for (CortexVertex cv : cg.vertexSet()) {
                    usedRois.put(cv.getCk(), true);
                }
            }

            pm.update();
        }

        return candidates;
    }

    @NotNull
    private Map<String, String> createContigEncoding() {
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
        return contigEncoding;
    }

    private String getContigCounts(KmerLookup kl, DepthFirstIterator<CortexVertex, CortexEdge> dfs, Map<CortexKmer, Boolean> usedRois, Map<String, String> contigEncoding) {
        StringBuilder sb = new StringBuilder();

        while (dfs.hasNext()) {
            CortexVertex cv = dfs.next();
            Set<Interval> loci = kl.findKmer(cv.getSk());
            if (usedRois.containsKey(cv.getCk())) {
                sb.append(".");
            } else if (loci.size() == 1) {
                for (Interval locus : loci) {
                    sb.append(contigEncoding.get(locus.getContig()));
                }
            } else {
                sb.append("_");
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
