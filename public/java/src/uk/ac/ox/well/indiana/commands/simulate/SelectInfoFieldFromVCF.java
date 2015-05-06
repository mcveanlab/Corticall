package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTable;

import java.io.PrintStream;
import java.util.HashSet;
import java.util.Map;

public class SelectInfoFieldFromVCF extends Module {
    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Argument(fullName="field", shortName="f", doc="Fields")
    public HashSet<String> FIELDS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        DataTable dt = new DataTable("attrs", "attrs");

        int index = 0;
        for (VariantContext vc : VCF) {
            Map<String, Object> attrs = vc.getAttributes();

            for (String field : FIELDS) {
                dt.set(String.valueOf(index), field, attrs.get(field));
            }

            index++;
        }

        out.println(dt);
    }
}
