package uk.ac.ox.well.cortexjdk.commands.discover.call;

import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;

import java.io.PrintStream;

public class FilterCalls extends Module {
    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VARIANTS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {

    }
}
