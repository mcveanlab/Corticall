package uk.ac.ox.well.indiana.attic.analyses.nahr;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;

public class SimulateVariation extends Module {
    @Argument(fullName="ref", shortName="r", doc="Reference genome")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="vcf", shortName="v", doc="Variants")
    public VCFFileReader VCF;

    @Argument(fullName="parent1", shortName="p1", doc="Parent 1")
    public String PARENT1;

    @Argument(fullName="parent2", shortName="p2", doc="Parent 2")
    public String PARENT2;

    @Output
    public PrintStream out;

    @Override
    public void execute() {


        log.info("Processing variants...");
        int numVariants = 0;
        for (VariantContext vc : VCF) {
            if (numVariants % 10000 == 0) {
                log.info("  {}", numVariants);
            }
            numVariants++;

            if (!vc.isFiltered()) {
                log.info("vc: {}", vc);


            }
        }
    }
}
