package uk.ac.ox.well.cortexjdk.commands.discover.eval;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

public class ViableKmerSizes extends Module {
    @Argument(fullName="fastq", shortName="f", doc="Fastqs")
    public ArrayList<File> FASTQS;

    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Argument(fullName="reference", shortName="R", doc="Reference")
    public IndexedFastaSequenceFile REF;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
    }
}
