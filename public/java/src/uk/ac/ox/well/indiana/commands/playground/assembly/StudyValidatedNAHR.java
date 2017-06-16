package uk.ac.ox.well.indiana.commands.playground.assembly;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by kiran on 16/06/2017.
 */
public class StudyValidatedNAHR extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROIs")
    public CortexGraph ROI;

    @Argument(fullName="lookup", shortName="l", doc="Lookup")
    public KmerLookup LOOKUP;

    @Argument(fullName="sequence", shortName="s", doc="Sequence")
    public FastaSequenceFile SEQUENCE;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<Integer> colors = new HashSet<>();
        colors.add(GRAPH.getColorForSampleName(CHILD));
        colors.addAll(GRAPH.getColorsForSampleNames(PARENTS));

        String seq = SEQUENCE.nextSequence().getBaseString();

        for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
            String sk = seq.substring(i, i + GRAPH.getKmerSize());
            CortexKmer ck = new CortexKmer(sk);
            CortexRecord cr = GRAPH.findRecord(ck);

            log.info("{}", recordToString(sk, cr, colors));
        }
    }

    private String recordToString(String sk, CortexRecord cr, Set<Integer> colors) {
        String kmer = cr.getKmerAsString();
        String cov = "";
        String ed = "";

        boolean fw = sk.equals(kmer);

        if (!fw) {
            kmer = SequenceUtils.reverseComplement(kmer);
        }

        int color = 0;
        for (int coverage : cr.getCoverages()) {
            if (colors.contains(color)) {
                cov += " " + coverage;
            }
            color++;
        }

        color = 0;
        for (String edge : cr.getEdgeAsStrings()) {
            if (colors.contains(color)) {
                ed += " " + (fw ? edge : SequenceUtils.reverseComplement(edge));
            }
            color++;
        }

        Set<String> lss = new TreeSet<>();
        if (LOOKUP != null) {
            Set<Interval> loci = LOOKUP.findKmer(kmer);

            if (loci != null && loci.size() > 0) {
                for (Interval locus : loci) {
                    String ls = locus.getContig() + ":" + locus.getStart() + "-" + locus.getEnd() + ":" + (locus.isPositiveStrand() ? "+" : "-");
                    lss.add(ls);
                }
            }
        }
        String lssCombined = Joiner.on(";").join(lss);

        return kmer + " " + cov + " " + ed + " " + lssCombined;
    }
}
