package uk.ac.ox.well.cortexjdk.commands.prefilter;

import htsjdk.samtools.SAMRecord;
import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Description;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.*;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContaminantStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.io.File;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

@Description(text="Find chains of contaminating kmers")
public class FindContamination extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="contamination", shortName="contam", doc="Contam")
    public CortexGraph CONTAM;

    @Argument(fullName = "drafts", shortName = "d", doc = "Drafts")
    public HashMap<String, IndexedReference> LOOKUPS;

    @Output
    public File out;

    //@Output(fullName="excluded_out", shortName="xo", doc="Excluded kmers output file")
    //public File contam_out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        Set<Integer> parentColors = new HashSet<>(GRAPH.getColorsForSampleNames(PARENTS));

        log.info("Colors:");
        log.info(" -   child: {}", GRAPH.getColorForSampleName(CHILD));
        log.info(" - parents: {}", GRAPH.getColorsForSampleNames(PARENTS));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Finding contamination")
                .message("records processed")
                .maxRecord(CONTAM.getNumRecords())
                .make(log);

        Map<CanonicalKmer, Boolean> roiKmers = new HashMap<>();
        for (CortexRecord rc : ROI) {
            roiKmers.put(rc.getCanonicalKmer(), false);
        }

        Set<CanonicalKmer> contamKmers = new HashSet<>();
        int numContamChains = 0;

        TraversalEngine e = new TraversalEngineFactory()
                .traversalDirection(BOTH)
                .combinationOperator(OR)
                .traversalColors(childColor)
                .joiningColors(parentColors)
                .stoppingRule(ContaminantStopper.class)
                .rois(ROI)
                .graph(GRAPH)
                .make();

        for (CortexRecord cr : CONTAM) {
            if (roiKmers.containsKey(cr.getCanonicalKmer()) && !roiKmers.get(cr.getCanonicalKmer())) {
                List<CortexVertex> l = e.walk(cr.getKmerAsString());

                List<String> pieces = new ArrayList<>();
                List<String> piece = new ArrayList<>();

                for (CortexVertex cv : l) {
                    String sk = cv.getKmerAsString();
                    CanonicalKmer ck = cv.getCanonicalKmer();

                    if (roiKmers.containsKey(ck)) {
                        if (piece.size() > 0) {
                            pieces.add(combineKmers(piece));
                            piece = new ArrayList<>();
                        }
                        roiKmers.put(ck, true);
                    } else {
                        piece.add(sk);
                    }
                }

                if (piece.size() > 0) {
                    pieces.add(combineKmers(piece));
                }


                boolean wellAligned = false;

                for (String p : pieces) {
                    for (String background : LOOKUPS.keySet()) {
                        List<SAMRecord> srs = LOOKUPS.get(background).getAligner().align(p);

                        int numAlignments = 0;
                        for (SAMRecord sr : srs) {
                            if (sr.getMappingQuality() > 0) {
                                numAlignments++;
                            }
                        }

                        if (numAlignments == 1) {
                            wellAligned = true;
                        }
                    }
                }

                for (CortexVertex v : l) {
                    if (roiKmers.containsKey(v.getCanonicalKmer())) {
                        roiKmers.put(v.getCanonicalKmer(), true);

                        if (!wellAligned) {
                            contamKmers.add(v.getCanonicalKmer());
                        }
                    }
                }

                if (wellAligned) {
                    numContamChains++;
                }
            }

            pm.update();
        }

        log.info("Found {} contamination kmer chains ({} kmers total)", numContamChains, contamKmers.size());

        log.info("Writing...");

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(makeCortexHeader());

        //CortexGraphWriter cgc = new CortexGraphWriter(contam_out);
        //cgc.setHeader(ROI.getHeader());

        int numKept = 0, numExcluded = 0;
        for (CortexRecord rr : ROI) {
            if (!contamKmers.contains(rr.getCanonicalKmer())) {
                //cgw.addRecord(rr);
                numKept++;
            } else {
                //cgc.addRecord(rr);
                cgw.addRecord(rr);
                numExcluded++;
            }
        }

        cgw.close();
        //cgc.close();

        log.info("  {}/{} ({}%) kept, {}/{} ({}%) excluded",
                numKept,     ROI.getNumRecords(), 100.0f * (float) numKept / (float) ROI.getNumRecords(),
                numExcluded, ROI.getNumRecords(), 100.0f * (float) numExcluded / (float) ROI.getNumRecords()
        );
    }

    @NotNull
    private CortexHeader makeCortexHeader() {
        CortexHeader ch = new CortexHeader();
        ch.setVersion(6);
        ch.setNumColors(1);
        ch.setKmerSize(GRAPH.getKmerSize());
        ch.setKmerBits(GRAPH.getKmerBits());

        CortexColor cc = new CortexColor();
        cc.setCleanedAgainstGraph(false);
        cc.setCleanedAgainstGraphName("");
        cc.setErrorRate(0);
        cc.setLowCovgKmersRemoved(false);
        cc.setLowCovgSupernodesRemoved(false);
        cc.setTipClippingApplied(false);
        cc.setTotalSequence(0);
        cc.setSampleName(CHILD);

        ch.addColor(cc);

        return ch;
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
