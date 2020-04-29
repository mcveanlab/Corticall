package uk.ac.ox.well.cortexjdk.playground.eval;

import com.google.common.base.Joiner;
import com.google.common.collect.Sets;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

public class FilterRefBasedVariants extends Module {
    @Argument(fullName="reference", shortName="R", doc="Reference")
    public ArrayList<FastaSequenceFile> REFERENCE;

    @Argument(fullName="childVariants", shortName="cv", doc="Child VCFs")
    public ArrayList<VCFFileReader> CHILD_VCFS;

    @Argument(fullName="parentVariants", shortName="pv", doc="Parent VCFs")
    public ArrayList<VCFFileReader> PARENT_VCFS;

    @Argument(fullName="window", shortName="w", doc="Window")
    public Integer WINDOW = 100;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        IntervalTreeMap<VariantContext> ivc = new IntervalTreeMap<>();
        for (VCFFileReader vcf : PARENT_VCFS) {
            for (VariantContext vc : vcf) {
                Interval interval = new Interval(vc.getContig(), vc.getStart(), vc.getEnd());

                ivc.put(interval, vc);
            }
        }

        for (VCFFileReader vcf : CHILD_VCFS) {
            for (VariantContext vc : vcf) {
                Interval interval = new Interval(vc.getContig(), vc.getStart() - WINDOW, vc.getEnd() + WINDOW);

                if (!vc.isFiltered() && !ivc.containsOverlapping(interval) && vc.getGenotype(0).isHomVar()) {
                    log.info("{}", vc);
                } else {
                    log.info("{}", vc);
                }
            }
        }
    }
}

