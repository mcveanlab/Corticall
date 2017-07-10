package uk.ac.ox.well.indiana.commands.playground.assemblies;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by kiran on 02/07/2017.
 */
public class CreateInheritanceTrack extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="lookups", shortName="l", doc="Lookup tables")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(new ArrayList<>(LOOKUPS.keySet()));

        for (String parent : LOOKUPS.keySet()) {
            int otherParentColor = -1;
            for (String p : LOOKUPS.keySet()) {
                if (!p.equals(parent)) {
                    otherParentColor = GRAPH.getColorForSampleName(p);
                }
            }

            ReferenceSequence rseq;
            while ((rseq = LOOKUPS.get(parent).getReferenceSequence().nextSequence()) != null) {
                String seq = rseq.getBaseString();

                int start = -1;
                int stop = -1;

                for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                    String sk = seq.substring(i, i + GRAPH.getKmerSize());
                    CortexKmer ck = new CortexKmer(sk);
                    CortexRecord cr = GRAPH.findRecord(ck);

                    boolean childHasCoverage = (cr != null && cr.getCoverage(childColor) > 0);
                    boolean otherParentHasNoCoverage = (cr != null && cr.getCoverage(otherParentColor) > 0);

                    if (childHasCoverage && !otherParentHasNoCoverage) {
                        if (start == -1) { start = i; }
                        stop = i + GRAPH.getKmerSize();
                    } else {
                        if (start != -1) {
                            log.info("{} {}:{}-{} {}", parent, rseq.getName().split("\\s+")[0], start, stop, stop - start);
                            start = -1;
                            stop = -1;
                        }
                    }
                }
            }
        }
    }
}
