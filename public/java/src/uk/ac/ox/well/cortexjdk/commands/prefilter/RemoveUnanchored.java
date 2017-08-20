package uk.ac.ox.well.cortexjdk.commands.prefilter;

import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.stoppingconditions.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.io.File;
import java.util.*;

/**
 * Created by kiran on 21/07/2017.
 */
public class RemoveUnanchored extends Module {
    @Argument(fullName = "graph", shortName = "g", doc = "Graph")
    public CortexGraph GRAPH;

    @Argument(fullName = "parents", shortName = "p", doc = "Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName = "child", shortName = "c", doc = "Child")
    public String CHILD;

    @Argument(fullName = "roi", shortName = "r", doc = "ROI")
    public CortexGraph ROI;

    @Argument(fullName = "drafts", shortName = "d", doc = "Drafts")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Output
    public File out;

    @Output(fullName = "unanchored_out", shortName = "uo", doc = "Unanchored output file")
    public File unanchored_out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        Set<Integer> parentColors = new HashSet<>(GRAPH.getColorsForSampleNames(PARENTS));
        //Set<Integer> ignoreColors = new HashSet<>(GRAPH.getColorsForSampleNames(IGNORE));

        log.info("Colors:");
        log.info(" -   child: {}", GRAPH.getColorForSampleName(CHILD));
        log.info(" - parents: {}", GRAPH.getColorsForSampleNames(PARENTS));
        //log.info(" -  ignore: {}", GRAPH.getColorsForSampleNames(IGNORE));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Finding unanchored kmers")
                .message("records processed")
                .maxRecord(ROI.getNumRecords())
                .make(log);

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(childColor)
                .joiningColors(parentColors)
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.OR)
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
                .rois(ROI)
                .stoppingRule(ContigStopper.class)
                .graph(GRAPH)
                .make();

        Set<CortexKmer> unanchored = new HashSet<>();
        int numUnanchoredChains = 0;

        Set<CortexKmer> rois = new HashSet<>();
        for (CortexRecord rr : ROI) {
            rois.add(rr.getCortexKmer());
        }

        for (CortexKmer rk : rois) {
            if (!unanchored.contains(rk)) {
                String contig = TraversalEngine.toContig(e.walk(rk.getKmerAsString()));
                StringBuilder annb = new StringBuilder();

                Set<CortexKmer> seenRois = new HashSet<>();

                for (int i = 0; i <= contig.length() - GRAPH.getKmerSize(); i++) {
                    String sk = contig.substring(i, i + GRAPH.getKmerSize());
                    CortexKmer ck = new CortexKmer(sk);

                    if (rois.contains(ck)) {
                        annb.append(".");

                        seenRois.add(ck);
                    } else {
                        String code = "?";
                        for (String background : LOOKUPS.keySet()) {
                            Set<Interval> intervals = LOOKUPS.get(background).findKmer(sk);

                            if (intervals.size() == 1) {
                                code = "1";
                            } else if (intervals.size() == 0 && !code.equals("1")) {
                                code = "?";
                            } else if (intervals.size() > 1 && !code.equals("1")) {
                                code = "_";
                            }
                        }

                        annb.append(code);
                    }
                }

                String ann = annb.toString();

                List<String> pieces = SequenceUtils.splitAtPositions(ann, SequenceUtils.computeSplits(ann, '.'));
                boolean leftAnchored = false;
                boolean novelsSeen = false;
                boolean rightAnchored = false;

                for (String piece : pieces) {
                    if (!piece.contains(".") && piece.contains("1")) {
                        if (!novelsSeen) {
                            leftAnchored = true;
                        } else {
                            rightAnchored = true;
                        }
                    }

                    if (piece.contains(".")) {
                        novelsSeen = true;
                    }
                }

                if (!leftAnchored || !rightAnchored) {
                    log.debug("    discard: {}", ann);

                    unanchored.addAll(seenRois);
                    numUnanchoredChains++;
                }
            }

            pm.update();
        }

        log.info("Found {} unanchored kmer chains ({} kmers total)", numUnanchoredChains, unanchored.size());

        log.info("Writing...");

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(ROI.getHeader());

        CortexGraphWriter cgt = new CortexGraphWriter(unanchored_out);
        cgt.setHeader(ROI.getHeader());

        int numKept = 0, numExcluded = 0;
        for (CortexRecord rr : ROI) {
            if (!unanchored.contains(rr.getCortexKmer())) {
                cgw.addRecord(rr);
                numKept++;
            } else {
                cgt.addRecord(rr);
                numExcluded++;
            }
        }

        cgw.close();
        cgt.close();

        log.info("  {}/{} ({}%) kept, {}/{} ({}%) excluded",
                numKept,     ROI.getNumRecords(), 100.0f * (float) numKept / (float) ROI.getNumRecords(),
                numExcluded, ROI.getNumRecords(), 100.0f * (float) numExcluded / (float) ROI.getNumRecords()
        );
    }
}
