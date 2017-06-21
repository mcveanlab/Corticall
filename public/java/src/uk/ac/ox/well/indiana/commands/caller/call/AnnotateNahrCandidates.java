package uk.ac.ox.well.indiana.commands.caller.call;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.util.*;

/**
 * Created by kiran on 21/06/2017.
 */
public class AnnotateNahrCandidates extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="sequences", shortName="s", doc="Sequences")
    public FastaSequenceFile SEQUENCES;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="refs", shortName="R", doc="References")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Argument(fullName="nahr", shortName="n", doc="NAHR")
    public FastaSequenceFile NAHR;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        List<Integer> parentColors = GRAPH.getColorsForSampleNames(PARENTS);
        List<Integer> recruitColors = GRAPH.getColorsForSampleNames(new ArrayList<>(LOOKUPS.keySet()));

        Set<Integer> colors = new TreeSet<>();
        colors.add(childColor);
        colors.addAll(parentColors);

        Map<CortexKmer, Set<String>> contigAssignments = new HashMap<>();

        List<String> backgrounds = new ArrayList<>(LOOKUPS.keySet());

        ReferenceSequence rseq;
        while ((rseq = SEQUENCES.nextSequence()) != null) {
            String seq = rseq.getBaseString();

            //log.info("{}", rseq.getName());

            String[] rname = rseq.getName().split("\\s+");
            String name = rname[0];

            for (int i = 0; i <= seq.length() - ROI.getKmerSize(); i++) {
                String sk = seq.substring(i, i + ROI.getKmerSize());

                ContainerUtils.add(contigAssignments, new CortexKmer(sk), name);

                /*
                List<String> intervals = new ArrayList<>();

                for (String background : backgrounds) {
                    Set<Interval> its = LOOKUPS.get(background).findKmer(sk);
                    if (its.size() == 1) {
                        Interval it = its.iterator().next();
                        intervals.add(background + ":" + it.getContig() + ":" + it.getStart() + "-" + it.getEnd() + ":" + (it.isPositiveStrand() ? "+" : "-"));
                    } else if (its.size() > 1) {
                        intervals.add(background + ":many");
                    } else {
                        intervals.add(background + ":none");
                    }
                }

                CortexRecord rr = ROI.findRecord(new CortexKmer(sk));
                CortexRecord cr = GRAPH.findRecord(new CortexKmer(sk));

                log.info("  {} {} {} {} {}", i, sk, rr != null, Joiner.on("\t").join(intervals), recordToString(sk, cr, colors));
                */
            }
        }

        log.info("NAHR:");

        while ((rseq = NAHR.nextSequence()) != null) {
            String seq = rseq.getBaseString();

            for (int i = 0; i <= seq.length() - ROI.getKmerSize(); i++) {
                String sk = seq.substring(i, i + ROI.getKmerSize());
                CortexKmer ck = new CortexKmer(sk);
                CortexRecord cr = GRAPH.findRecord(ck);

                log.info("  {} {} {}", sk, contigAssignments.containsKey(ck) ? contigAssignments.get(ck) : null, recordToString(sk, cr, colors));

            }
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

        /*
        Set<String> lss = new TreeSet<>();
        if (LOOKUPS != null) {
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
        */

        return kmer + " " + cov + " " + ed;
    }
}
