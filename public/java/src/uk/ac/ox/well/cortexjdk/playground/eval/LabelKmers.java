package uk.ac.ox.well.cortexjdk.playground.eval;

import com.google.common.base.Joiner;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.containers.ContainerUtils;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;

import java.io.PrintStream;
import java.util.*;

public class LabelKmers extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="rois", shortName="r", doc="ROIs")
    public CortexGraph ROIS;

    @Argument(fullName="homVarKmers", shortName="v", doc="Homvar kmers")
    public CortexGraph HOM_VAR_KMERS;

    @Argument(fullName="backgrounds", shortName="b", doc="Background labels")
    public LinkedHashSet<String> BACKGROUNDS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CanonicalKmer, Boolean> kmerGenotypes = new HashMap<>();
        Map<CanonicalKmer, String> kmerLabels = new HashMap<>();

        for (CortexRecord cr : ROIS) {
            kmerGenotypes.put(cr.getCanonicalKmer(), false);

            List<String> labels = new ArrayList<>();

            for (String bg : BACKGROUNDS) {
                CortexRecord gr = GRAPH.findRecord(cr.getCanonicalKmer());
                if (gr != null && gr.getCoverage(GRAPH.getColorForSampleName(bg)) > 0) {
                    labels.add(bg);
                }
            }

            if (labels.size() > 0) {
                kmerLabels.put(cr.getCanonicalKmer(), labels.get(0));
            } else {
                kmerLabels.put(cr.getCanonicalKmer(), "unknown");
            }
        }

        for (CortexRecord cr : HOM_VAR_KMERS) {
            kmerGenotypes.put(cr.getCanonicalKmer(), true);
        }

        Map<String, Integer> countRef = new HashMap<>();
        Map<String, Integer> countVar = new HashMap<>();

        for (String bg : BACKGROUNDS) {
            countRef.put(bg, 0);
            countVar.put(bg, 0);
        }

        for (CanonicalKmer ck : kmerGenotypes.keySet()) {
            String label = kmerLabels.get(ck);
            boolean genotype = kmerGenotypes.get(ck);

            if (!genotype) { // (hom-ref)
                ContainerUtils.increment(countRef, label);
            } else {
                ContainerUtils.increment(countVar, label);
            }
        }

        for (String bg : BACKGROUNDS) {
            out.println(Joiner.on("\t").join(bg, countRef.get(bg), countVar.get(bg)));
        }
    }
}

