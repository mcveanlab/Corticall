package uk.ac.ox.well.indiana.commands.gg;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTable;
import uk.ac.ox.well.indiana.utils.containers.DataTables;

import java.io.File;
import java.io.PrintStream;

public class CompareGraphCallsToVCF extends Module {
    @Argument(fullName="graphCalls", shortName="gc", doc="Graph calls")
    public File GRAPH_CALLS;

    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Argument(fullName="ref", shortName="r", doc="Reference")
    public IndexedFastaSequenceFile REF;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        DataTables dts = new DataTables(GRAPH_CALLS);
        DataTable dt = dts.getTable("variantCalls");

        for (VariantContext vc : VCF) {
            log.info("vc: {}", vc);

            String chr = vc.getChr();
            int start = vc.getStart();
            int end = vc.getEnd();
        }
    }
}
