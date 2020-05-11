package uk.ac.ox.well.cortexjdk.playground.eval;

import com.google.common.base.Joiner;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

import static htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType.VCF;

public class FilterRefBasedVariants extends Module {
    @Argument(fullName="reference", shortName="R", doc="Reference")
    public ArrayList<FastaSequenceFile> REFERENCE;

    @Argument(fullName="childVariants", shortName="cv", doc="Child VCFs")
    public ArrayList<VCFFileReader> CHILD_VCFS;

    @Argument(fullName="parentVariants", shortName="pv", doc="Parent VCFs", required=false)
    public ArrayList<VCFFileReader> PARENT_VCFS;

    @Argument(fullName="window", shortName="w", doc="Window")
    public Integer WINDOW = 10000;

    @Output
    public File out;

    @Override
    public void execute() {
        IntervalTreeMap<VariantContext> ivc = new IntervalTreeMap<>();
        if (PARENT_VCFS != null) {
            for (VCFFileReader vcf : PARENT_VCFS) {
                for (VariantContext vc : vcf) {
                    if (!vc.isFiltered() && (!vc.hasGenotypes() || !vc.getGenotype(0).isHomRef())) {
                        Interval interval = new Interval(vc.getContig(), vc.getStart(), vc.getEnd());

                        if (interval.length() < 2*WINDOW) {
                            ivc.put(interval, vc);
                        }
                    }
                }
            }
        }

        VariantContextWriter vcw = new VariantContextWriterBuilder()
                .setOutputFile(out)
                .setOutputFileType(VCF)
                .setOption(Options.DO_NOT_WRITE_GENOTYPES)
                .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .unsetOption(Options.INDEX_ON_THE_FLY)
                .build();

        VCFHeader vcfHeader = new VCFHeader();
        SAMSequenceDictionary sd = buildMergedSequenceDictionary(REFERENCE);
        vcfHeader.setSequenceDictionary(sd);
        vcw.writeHeader(vcfHeader);

        for (VCFFileReader vcf : CHILD_VCFS) {
            for (VariantContext vc : vcf) {
                Interval interval = new Interval(vc.getContig(), vc.getStart() - WINDOW, vc.getEnd() + WINDOW);

                if (!vc.isFiltered() && (!vc.hasGenotypes() || !vc.getGenotype(0).isHomRef()) && !ivc.containsOverlapping(interval)) {
                    log.info("{}", vc);

                    vcw.add(vc);
                }
            }
        }

        vcw.close();
    }

    @NotNull
    private SAMSequenceDictionary buildMergedSequenceDictionary(List<FastaSequenceFile> refs) {
        List<SAMSequenceRecord> ssrs = new ArrayList<>();

        for (FastaSequenceFile ref : refs) {
            ssrs.addAll(ref.getSequenceDictionary().getSequences());
        }

        return new SAMSequenceDictionary(ssrs);
    }
}

