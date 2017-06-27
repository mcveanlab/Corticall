package uk.ac.ox.well.indiana.commands.caller.call;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.jetbrains.annotations.NotNull;
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
 * Created by kiran on 27/06/2017.
 */
public class AnnotateContigs extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="drafts", shortName="d", doc="Drafts")
    public TreeMap<String, KmerLookup> LOOKUPS;

    @Argument(fullName="sequences", shortName="s", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(PARENTS);
        List<Integer> recruitColors = GRAPH.getColorsForSampleNames(new ArrayList<>(LOOKUPS.keySet()));
        int refColor = GRAPH.getColorForSampleName("ref");

        out.println(Joiner.on("\t").join("name", "sk", "ck", "cov_" + CHILD, "cov_" + Joiner.on("\tcov_").join(PARENTS), "cov_" + Joiner.on("\tcov_").join(LOOKUPS.keySet()), "cov_ref", "is_novel", "is_filled_gap", "is_recovered_kmer", LOOKUPS.keySet()));

        List<ReferenceSequence> rseqs = new ArrayList<>();
        ReferenceSequence aseq;
        while ((aseq = CONTIGS.nextSequence()) != null) {
            rseqs.add(aseq);
        }

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Annotating contigs")
                .maxRecord(rseqs.size())
                .message("contigs annotated")
                .make(log);

        for (ReferenceSequence rseq : rseqs) {
            String[] pieces = rseq.getName().split("\\s+");
            String seq = rseq.getBaseString();

            for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                String sk = seq.substring(i, i + GRAPH.getKmerSize());
                CortexKmer ck = new CortexKmer(sk);
                CortexRecord cr = GRAPH.findRecord(ck);

                List<Integer> coverages = new ArrayList<>();
                coverages.add(cr != null ? cr.getCoverage(childColor) : 0);
                parentColors.forEach(c -> { coverages.add(cr != null ? cr.getCoverage(c) : 0); });
                recruitColors.forEach(c -> { coverages.add(cr != null ? cr.getCoverage(c) : 0); });
                coverages.add(cr != null ? cr.getCoverage(refColor) : 0);

                List<String> allIntervals = new ArrayList<>();
                for (String background : LOOKUPS.keySet()) {
                    Set<Interval> its = new TreeSet<>(LOOKUPS.get(background).findKmer(sk));
                    List<String> intervalStrings = new ArrayList<>();

                    if (its.size() > 0) {
                        for (Interval it : its) {
                            intervalStrings.add(it.getContig() + ":" + it.getStart() + "-" + it.getEnd() + ":" + (it.isPositiveStrand() ? "+" : "-"));
                        }
                    } else {
                        intervalStrings.add("NA");
                    }

                    allIntervals.add(Joiner.on(";").join(intervalStrings));
                }

                CortexRecord rr = ROI.findRecord(ck);

                boolean isNovel = rr != null;
                boolean isFilledGap = cr != null && cr.getCoverage(childColor) == 0;
                boolean isRecovered = cr != null && cr.getCoverage(childColor) < GRAPH.getColor(childColor).getLowCovSupernodesThreshold();

                out.println(Joiner.on("\t").join(pieces[0], sk, ck, Joiner.on("\t").join(coverages), isNovel, isFilledGap, isRecovered, Joiner.on("\t").join(allIntervals)));
            }

            pm.update();
        }
    }
}
