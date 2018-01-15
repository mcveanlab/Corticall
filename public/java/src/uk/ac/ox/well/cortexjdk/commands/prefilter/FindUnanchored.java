package uk.ac.ox.well.cortexjdk.commands.prefilter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Description;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraphWriter;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.File;
import java.util.*;

/**
 * Created by kiran on 21/07/2017.
 */
@Description(text="Find novel kmers in contigs that cannot be confidently placed on any available reference")
public class FindUnanchored extends Module {
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

    //@Output(fullName = "excluded_out", shortName = "xo", doc = "Excluded kmers output file")
    //public File unanchored_out;

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

        Set<CanonicalKmer> unanchored = new HashSet<>();
        int numUnanchoredChains = 0;

        Set<CanonicalKmer> rois = new HashSet<>();
        for (CortexRecord rr : ROI) {
            rois.add(rr.getCanonicalKmer());
        }

        for (CanonicalKmer rk : rois) {
            if (!unanchored.contains(rk)) {
                String contig = TraversalEngine.toContig(e.walk(rk.getKmerAsString()));
                StringBuilder annb = new StringBuilder();

                Set<CanonicalKmer> seenRois = new HashSet<>();

                for (int i = 0; i <= contig.length() - GRAPH.getKmerSize(); i++) {
                    String sk = contig.substring(i, i + GRAPH.getKmerSize());
                    CanonicalKmer ck = new CanonicalKmer(sk);

                    if (rois.contains(ck)) {
                        annb.append(".");

                        seenRois.add(ck);
                    } else {
                        annb.append("_");
                    }
                }

                String ann = annb.toString();

                List<String> pieces = SequenceUtils.splitAtPositions(ann, SequenceUtils.computeSplits(ann, '.'));
                boolean novelsSeen = false;
                boolean hasAlignments = false;

                for (int i = 0; i < pieces.size(); i++) {
                    if (pieces.contains(".")) {
                        novelsSeen = true;
                    } else {
                        for (String background : LOOKUPS.keySet()) {
                            List<SAMRecord> srs = LOOKUPS.get(background).getAligner().align(pieces.get(i));

                            for (SAMRecord sr : srs) {
                                if (sr.getMappingQuality() > 0) {
                                    hasAlignments = true;
                                    break;
                                }
                            }
                        }
                    }

                    if (novelsSeen && hasAlignments) {
                        break;
                    }
                }

                if (!hasAlignments) {
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

        int numKept = 0, numExcluded = 0;
        for (CortexRecord rr : ROI) {
            if (!unanchored.contains(rr.getCanonicalKmer())) {
                numKept++;
            } else {
                cgw.addRecord(rr);
                numExcluded++;
            }
        }

        cgw.close();

        log.info("  {}/{} ({}%) kept, {}/{} ({}%) excluded",
                numKept,     ROI.getNumRecords(), 100.0f * (float) numKept / (float) ROI.getNumRecords(),
                numExcluded, ROI.getNumRecords(), 100.0f * (float) numExcluded / (float) ROI.getNumRecords()
        );
    }
}
