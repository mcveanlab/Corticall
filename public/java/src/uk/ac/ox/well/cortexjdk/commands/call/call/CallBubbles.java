package uk.ac.ox.well.cortexjdk.commands.call.call;

import com.google.common.base.Joiner;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.caller.Bubble;
import uk.ac.ox.well.cortexjdk.utils.caller.BubbleCaller;
import uk.ac.ox.well.cortexjdk.utils.caller.BubbleCallerFactory;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraphWriter;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by kiran on 29/08/2017.
 */
public class CallBubbles extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="links", shortName="l", doc="Links")
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="references", shortName="R", doc="References")
    public HashMap<String, KmerLookup> REFERENCES;

    @Output
    public PrintStream out;

    @Output(fullName="rout", shortName="ro", doc="Remaining ROI out")
    public File rout;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(new ArrayList<>(REFERENCES.keySet()));

        log.info("Colors:");
        log.info("  {}", childColor);
        log.info("  {}", parentColors);

        BubbleCaller bc = new BubbleCallerFactory()
                .alternateColor(childColor)
                .referenceColors(parentColors)
                .references(REFERENCES)
                .graph(GRAPH)
                .links(LINKS)
                .rois(ROI)
                .make();

        out.println(Joiner.on("\t").join("contig", "start", "type", "ref", "alt", "flank5p", "flank3p", "nkCount", "nk", "nks"));

        Map<CortexKmer, Boolean> used = loadRois();
        int numBubbles = 0;
        int numNovelKmersInVariants = 0;

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing novel kmers")
                .message("records processed")
                .maxRecord(used.size())
                .make(log);

        for (CortexKmer rk : used.keySet()) {
            if (!used.get(rk)) {
                Set<Bubble> bs = bc.call(rk.getKmerAsString());
            }

            pm.update();
        }

        log.info("Found {} bubbles.  Used {}/{} novel kmers", numBubbles, numNovelKmersInVariants, used.size());

        CortexGraphWriter cgw = new CortexGraphWriter(rout);
        cgw.setHeader(ROI.getHeader());

        for (CortexRecord rr : ROI) {
            if (!used.get(rr.getCortexKmer())) {
                cgw.addRecord(rr);
            }
        }

        cgw.close();
    }

    private Map<CortexKmer, Boolean> loadRois() {
        Map<CortexKmer, Boolean> rrs = new HashMap<>();
        for (CortexRecord rr : ROI) {
            rrs.put(rr.getCortexKmer(), false);
        }
        return rrs;
    }
}
