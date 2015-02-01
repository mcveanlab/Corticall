package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;
import java.util.ArrayList;

public class FindVariantsInContigs extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Argument(fullName="reference", shortName="r", doc="Reference")
    public ArrayList<IndexedFastaSequenceFile> REFS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
    int lastPos = 0;
        for (VariantContext vc : VCF) {
            String chr = vc.getChr();
            int start = vc.getStart();
            int end = vc.getEnd();

            log.info("{} {} {} {} {}", chr, start, end, start - lastPos, vc);

            lastPos = end;
        }

    }
}
