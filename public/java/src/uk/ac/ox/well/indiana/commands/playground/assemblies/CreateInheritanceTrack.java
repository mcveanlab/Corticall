package uk.ac.ox.well.indiana.commands.playground.assemblies;

import com.google.common.base.Joiner;
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

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="reference", shortName="r", doc="Reference")
    public KmerLookup REFERENCE;

    /*
    @Argument(fullName="lookups", shortName="l", doc="Lookup tables")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Argument(fullName="reference", shortName="r", doc="Reference genome")
    public FastaSequenceFile REF;

    @Argument(fullName="mappings", shortName="m", doc="Coordinate mapping tables")
    public HashMap<String, File> MAPPING;
    */

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(PARENTS);

        Map<Integer, Integer> colorToIndexMap = new HashMap<>();
        for (int i = 0; i < parentColors.size(); i++) {
            colorToIndexMap.put(parentColors.get(i), i);
        }
        //List<Integer> recruitColors = GRAPH.getColorsForSampleNames(new ArrayList<>(LOOKUPS.keySet()));

        Map<Interval, Integer> kmerInheritance = new TreeMap<>();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing records...")
                .maxRecord(GRAPH.getNumRecords())
                .message("records processed")
                .make(log);

        for (CortexRecord cr : GRAPH) {
            boolean childHasCoverage = cr.getCoverage(childColor) > 0;
            int numParentsWithCoverage = 0;
            int parentColorWithCoverage = -1;

            for (int c : parentColors) {
                if (cr.getCoverage(c) > 0) {
                    numParentsWithCoverage++;
                    parentColorWithCoverage = c;
                }
            }

            if (childHasCoverage && numParentsWithCoverage == 1) {
                Set<Interval> loci = REFERENCE.findKmer(cr.getKmerAsString());

                if (loci.size() == 1) {
                    kmerInheritance.put(loci.iterator().next(), parentColorWithCoverage);
                }
            }

            pm.update();
        }

        log.info("  {} useful records found", kmerInheritance.size());

        log.info("Writing...");
        for (Interval it : kmerInheritance.keySet()) {
            out.println(Joiner.on(" ").join(it.getContig(), it.getStart(), it.getEnd(), colorToIndexMap.get(kmerInheritance.get(it))));
        }
    }
}
