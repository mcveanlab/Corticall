package uk.ac.ox.well.cortexjdk.commands.inheritance;

import com.google.common.base.Joiner;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/**
 * Created by kiran on 09/08/2017.
 */
public class VCFToInheritanceTrack extends Module {
    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Argument(fullName="parent", shortName="p", doc="Parent")
    public TreeSet<String> PARENTS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        List<String> parentNamesInVCF = getParentNamesInVCF();
        Set<String> headerSampleNames = new TreeSet<>(VCF.getFileHeader().getGenotypeSamples());

        out.println(Joiner.on("\t").join("#CHROM", "POS", "REF", "ALT", Joiner.on("\t").join(headerSampleNames)));

        for (VariantContext vc : VCF) {
            Genotype g0 = vc.getGenotype(parentNamesInVCF.get(0));
            Genotype g1 = vc.getGenotype(parentNamesInVCF.get(1));

            List<Integer> ids = new ArrayList<>();
            for (String sn : headerSampleNames) {
                Genotype gs = vc.getGenotype(sn);

                if (g0.getGenotypeString().equals(gs.getGenotypeString()) && !g1.getGenotypeString().equals(gs.getGenotypeString())) {
                    ids.add(0);
                } else {
                    ids.add(1);
                }
            }

            out.println(Joiner.on("\t").join(vc.getContig(), vc.getStart(), vc.getReference().toString(), vc.getAltAlleleWithHighestAlleleCount().toString(), Joiner.on("\t").join(ids)));
        }
    }

    private List<String> getParentNamesInVCF() {
        List<String> parentNamesInVCF = new ArrayList<>();
        for (String sampleNameFromVCF : VCF.getFileHeader().getGenotypeSamples()) {
            for (String parentNameFromUser : PARENTS) {
                if (sampleNameFromVCF.contains(parentNameFromUser)) {
                    parentNamesInVCF.add(sampleNameFromVCF);
                }
            }
        }
        return parentNamesInVCF;
    }
}
