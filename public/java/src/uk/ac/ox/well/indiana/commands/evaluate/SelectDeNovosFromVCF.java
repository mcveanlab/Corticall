package uk.ac.ox.well.indiana.commands.evaluate;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;

public class SelectDeNovosFromVCF extends Module {
    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="mother", shortName="m", doc="Mother")
    public String MOTHER;

    @Argument(fullName="father", shortName="f", doc="Father")
    public String FATHER;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        for (VariantContext vc : VCF) {
            if ( vc.getGenotype(CHILD).isCalled() && vc.getGenotype(MOTHER).isCalled() && vc.getGenotype(FATHER).isCalled() &&
                !vc.getGenotype(CHILD).getType().equals(vc.getGenotype(MOTHER).getType()) &&
                !vc.getGenotype(CHILD).getType().equals(vc.getGenotype(FATHER).getType()) &&
                 vc.getGenotype(CHILD).getGQ() > 90 && vc.getGenotype(MOTHER).getGQ() > 90 && vc.getGenotype(FATHER).getGQ() > 90 &&
                 vc.getGenotype(CHILD).getDP() > 50 && vc.getGenotype(MOTHER).getDP() > 50 && vc.getGenotype(FATHER).getDP() > 50
                ) {
                out.println(vc);
            }
        }
    }
}
