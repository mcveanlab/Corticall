package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTable;

import java.io.PrintStream;
import java.util.ArrayList;

public class CountDeNovoVariantsInVCF extends Module {
    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Argument(fullName="parent", shortName="p", doc="Parent")
    public ArrayList<String> PARENTS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        String s0 = PARENTS.get(0);
        String s1 = PARENTS.get(1);

        DataTable dt = new DataTable("denovo", "denovo count");

        for (VariantContext vc : VCF) {
            Genotype g0 = vc.getGenotype(s0);
            Genotype g1 = vc.getGenotype(s1);

            Allele a0 = g0.getAllele(0);
            Allele a1 = g1.getAllele(0);

            for (Genotype g : vc.getGenotypes()) {
                if (!g.getSampleName().equals(s0) && !g.getSampleName().equals(s1)) {
                    Allele a = g.getAllele(0);

                    if (!a0.equals(a) && !a1.equals(a)) {
                        dt.increment("denovo", g.getSampleName());
                    }
                }
            }
        }

        out.println(dt);
    }
}
