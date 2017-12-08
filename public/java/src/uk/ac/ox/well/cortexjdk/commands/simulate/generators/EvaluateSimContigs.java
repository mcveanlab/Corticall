package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMRecord;
import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngine;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineFactory;

import java.util.*;

public class EvaluateSimContigs extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Links")
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="ref", shortName="R", doc="Ref")
    public KmerLookup REF;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(ROI.getColor(0).getSampleName());

        TraversalEngine er = new TraversalEngineFactory()
                .traversalColor(childColor)
                .graph(GRAPH)
                .make();

        TraversalEngine el = new TraversalEngineFactory()
                .traversalColor(childColor)
                .graph(GRAPH)
                .links(LINKS)
                .make();

        Map<CanonicalKmer, String> seenStrs = new HashMap<>();
        Map<CanonicalKmer, Boolean> seen = new HashMap<>();
        for (CortexRecord rr : ROI) {
            seenStrs.put(rr.getCanonicalKmer(), null);
            seen.put(rr.getCanonicalKmer(), false);
        }

        for (CortexRecord rr : ROI) {
            if (seenStrs.containsKey(rr.getCanonicalKmer()) && seenStrs.get(rr.getCanonicalKmer()) != null) {
                log.info("{} {}", rr.getCanonicalKmer(), seenStrs.get(rr.getCanonicalKmer()));
            } else {
                List<CortexVertex> erw = er.walk(rr.getKmerAsString());
                List<CortexVertex> elw = el.walk(rr.getKmerAsString());
                //List<CortexVertex> ell = longWalk(seen, el, rr.getCanonicalKmer());

                List<SAMRecord> srw = REF.getAligner().align(TraversalEngine.toContig(erw));
                List<SAMRecord> slw = REF.getAligner().align(TraversalEngine.toContig(elw));

                String out = Joiner.on(" ").join(rr.getKmerAsString(), erw.size(), elw.size());
                log.info("{} {}", rr.getCanonicalKmer(), out);

                log.info("  - srw: {}", srw);
                log.info("  - slw: {}", slw);

                seenStrs.put(rr.getCanonicalKmer(), out);
                for (CortexVertex cv : erw) {
                    if (seenStrs.containsKey(cv.getCanonicalKmer()) && seenStrs.get(cv.getCanonicalKmer()) == null) {
                        seenStrs.put(cv.getCanonicalKmer(), out);
                    }
                }
            }
        }
    }

    @NotNull
    private List<CortexVertex> longWalk(Map<CanonicalKmer, Boolean> seen, TraversalEngine e, CanonicalKmer ck) {
        List<CortexVertex> w = e.walk(ck.getKmerAsString());

        //int wOldSize = w.size();

        boolean extended = false;
        do {
            extended = false;
            List<List<CortexVertex>> extFwd = new ArrayList<>();

            Set<CortexVertex> nvs = e.getNextVertices(w.get(w.size() - 1).getKmerAsByteKmer());
            for (CortexVertex cv : nvs) {
                List<CortexVertex> wn = e.walk(cv.getKmerAsString(), true);
                wn.add(0, cv);

                boolean hasNovels = false;

                for (CortexVertex v : wn) {
                    if (seen.containsKey(v.getCanonicalKmer()) && !seen.get(v.getCanonicalKmer())) {
                        hasNovels = true;
                        break;
                    }
                }

                if (hasNovels) {
                    extFwd.add(wn);
                }
            }

            if (extFwd.size() == 1) {
                w.addAll(extFwd.get(0));
                extended = true;

                for (CortexVertex v : extFwd.get(0)) {
                    if (seen.containsKey(v.getCanonicalKmer())) {
                        seen.put(v.getCanonicalKmer(), true);
                    }
                }
            }
        } while (extended);

        do {
            extended = false;
            List<List<CortexVertex>> extRev = new ArrayList<>();

            Set<CortexVertex> pvs = e.getPrevVertices(w.get(0).getKmerAsByteKmer());
            for (CortexVertex cv : pvs) {
                List<CortexVertex> wp = e.walk(cv.getKmerAsString(), false);
                wp.add(cv);

                boolean hasNovels = false;

                for (CortexVertex v : wp) {
                    if (seen.containsKey(v.getCanonicalKmer()) && !seen.get(v.getCanonicalKmer())) {
                        hasNovels = true;
                        break;
                    }
                }

                if (hasNovels) {
                    extRev.add(wp);
                }
            }

            if (extRev.size() == 1) {
                w.addAll(0, extRev.get(0));
                extended = true;

                for (CortexVertex v : extRev.get(0)) {
                    if (seen.containsKey(v.getCanonicalKmer())) {
                        seen.put(v.getCanonicalKmer(), true);
                    }
                }
            }
        } while (extended);

        //int wNewSize = w.size();

        return w;
    }
}
