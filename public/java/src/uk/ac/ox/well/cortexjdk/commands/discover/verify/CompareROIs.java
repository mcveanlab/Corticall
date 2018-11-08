package uk.ac.ox.well.cortexjdk.commands.discover.verify;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class CompareROIs extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="truth", shortName="t", doc="Truth")
    public CortexGraph TRUTH;

    @Argument(fullName="eval", shortName="e", doc="Eval")
    public CortexGraph EVAL;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CanonicalKmer, CortexRecord> trs = new HashMap<>();
        Map<CanonicalKmer, CortexRecord> ers = new HashMap<>();
        Set<CanonicalKmer> all = new HashSet<>();

        for (CortexRecord tr : TRUTH) {
            trs.put(tr.getCanonicalKmer(), tr);
            all.add(tr.getCanonicalKmer());
        }

        for (CortexRecord er : EVAL) {
            ers.put(er.getCanonicalKmer(), er);
            all.add(er.getCanonicalKmer());
        }

        int privateToTruth = 0, privateToEval = 0, overlap = 0;

        Set<CanonicalKmer> kmersPrivateToEval = new HashSet<>();

        for (CanonicalKmer ck : all) {
            if (trs.containsKey(ck)) {
                if (ers.containsKey(ck)) {
                    overlap++;
                } else {
                    privateToTruth++;

                    log.debug("{} {} {}", EVAL.getSampleName(0), GRAPH.getColorForSampleName(EVAL.getSampleName(0)), GRAPH.findRecord(ck));
                }
            } else {
                privateToEval++;
                kmersPrivateToEval.add(ck);
            }
        }

        log.info("t={} e={} pt={} pe={} o={}", trs.size(), ers.size(), privateToTruth, privateToEval, overlap);
        out.println("t=" + trs.size() + " e=" + ers.size() + " pt=" + privateToTruth + " pe=" + privateToEval + " o=" + overlap);

        for (CanonicalKmer ck : kmersPrivateToEval) {
            log.info("pe: {} {} {}", ck, GRAPH.findRecord(ck), SequenceUtils.computeCompressionRatio(ck));
        }
    }
}
