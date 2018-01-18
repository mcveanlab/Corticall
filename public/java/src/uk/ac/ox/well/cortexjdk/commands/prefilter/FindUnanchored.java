package uk.ac.ox.well.cortexjdk.commands.prefilter;

import htsjdk.samtools.SAMRecord;
import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Description;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraphWriter;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
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

    @Argument(fullName = "links", shortName = "l", doc = "Links", required = false)
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName = "parents", shortName = "p", doc = "Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName = "child", shortName = "c", doc = "Child")
    public String CHILD;

    @Argument(fullName = "roi", shortName = "r", doc = "ROI")
    public CortexGraph ROI;

    @Argument(fullName = "drafts", shortName = "d", doc = "Drafts")
    public HashMap<String, IndexedReference> LOOKUPS;

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
                .links(LINKS)
                .make();

        Set<CanonicalKmer> used = new HashSet<>();
        Set<CanonicalKmer> unanchored = new HashSet<>();
        int numUnanchoredChains = 0;

        Set<CanonicalKmer> rois = new HashSet<>();
        for (CortexRecord rr : ROI) {
            rois.add(rr.getCanonicalKmer());
        }

        for (CanonicalKmer rk : rois) {
            if (!used.contains(rk)) {
                String contig = TraversalEngine.toContig(e.walk(rk.getKmerAsString()));

                Set<CanonicalKmer> seenRois = new HashSet<>();
                List<String> pieces = new ArrayList<>();
                List<String> piece = new ArrayList<>();

                for (int i = 0; i <= contig.length() - GRAPH.getKmerSize(); i++) {
                    String sk = contig.substring(i, i + GRAPH.getKmerSize());
                    CanonicalKmer ck = new CanonicalKmer(sk);

                    if (rois.contains(ck)) {
                        if (piece.size() > 0) {
                            pieces.add(combineKmers(piece));
                            piece = new ArrayList<>();
                        }
                        seenRois.add(ck);
                    } else {
                        piece.add(sk);
                    }
                }

                if (piece.size() > 0) {
                    pieces.add(combineKmers(piece));
                }

                boolean hasAlignments = false;
                for (String p : pieces) {
                    for (String background : LOOKUPS.keySet()) {
                        List<SAMRecord> srs = LOOKUPS.get(background).getAligner().align(p);

                        for (SAMRecord sr : srs) {
                            int numAlignments = 0;
                            if (sr.getMappingQuality() > 0) {
                                numAlignments++;
                            }

                            if (numAlignments == 1) {
                                hasAlignments = true;
                                break;
                            }
                        }
                    }

                    if (hasAlignments) {
                        break;
                    }
                }

                if (!hasAlignments) {
                    unanchored.addAll(seenRois);
                    numUnanchoredChains++;
                }
                used.addAll(seenRois);
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

    @NotNull
    private String combineKmers(List<String> piece) {
        StringBuilder sb = new StringBuilder();
        for (String s : piece) {
            if (sb.length() == 0) {
                sb.append(s);
            } else {
                sb.append(s.substring(s.length() - 1, s.length()));
            }
        }
        return sb.toString();
    }
}
