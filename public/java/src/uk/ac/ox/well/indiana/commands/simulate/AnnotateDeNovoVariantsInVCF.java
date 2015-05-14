package uk.ac.ox.well.indiana.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTable;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;

public class AnnotateDeNovoVariantsInVCF extends Module {
    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Argument(fullName="parent", shortName="p", doc="Parent")
    public ArrayList<String> PARENTS;

    @Output
    public File out;

    @Output(fullName="statsOut", shortName="so")
    public PrintStream sout;

    @Override
    public void execute() {
        String s0 = PARENTS.get(0);
        String s1 = PARENTS.get(1);

        for (String s : VCF.getFileHeader().getGenotypeSamples()) {
            if (s.contains(s0)) { s0 = s; }
            if (s.contains(s1)) { s1 = s; }
        }

        DataTable dt = new DataTable("denovo", "denovo count");

        VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder()
                .setOutputFile(out)
                .setReferenceDictionary(VCF.getFileHeader().getSequenceDictionary());

        VariantContextWriter vcw = vcwb.build();

        VCFHeader header = VCF.getFileHeader();
        header.addMetaDataLine(new VCFInfoHeaderLine("SAMPLES_WITH_DENOVOS", 1, VCFHeaderLineType.String, "Comma-separated list of samples with de novo variants"));
        header.addMetaDataLine(new VCFInfoHeaderLine("SIMID", 1, VCFHeaderLineType.String, "The ID of the simulated variant"));

        vcw.writeHeader(header);

        int vindex = 0;
        for (VariantContext vc : VCF) {
            if (!vc.isFiltered()) {
                Genotype g0 = vc.getGenotype(s0);
                Genotype g1 = vc.getGenotype(s1);

                Allele a0 = g0.getAllele(0);
                Allele a1 = g1.getAllele(0);

                boolean hasDeNovoVariants = false;
                Set<String> samplesWithDeNovoVariants = new TreeSet<String>();

                for (Genotype g : vc.getGenotypes()) {
                    if (!g.getSampleName().contains(s0) && !g.getSampleName().contains(s1)) {
                        Allele a = g.getAllele(0);

                        //if (a.isCalled() && a0.isCalled() && a1.isCalled() && a.isNonReference() && !a0.equals(a) && !a1.equals(a)) {
                        if (a.isCalled() && a0.isCalled() && a1.isCalled() && !a0.equals(a) && !a1.equals(a)) {
                            dt.set(g.getSampleName(), "sample", g.getSampleName());
                            dt.set(g.getSampleName(), "isFilteredOut", vc.isFiltered());
                            dt.increment(g.getSampleName(), vc.getType().name());

                            hasDeNovoVariants = true;
                            samplesWithDeNovoVariants.add(g.getSampleName());
                        }
                    }
                }

                if (hasDeNovoVariants) {
                    VariantContext newvc = (new VariantContextBuilder(vc))
                            .attribute("SAMPLES_WITH_DENOVOS", Joiner.on(",").join(samplesWithDeNovoVariants))
                            .make();

                    vc = newvc;
                }
            }

            VariantContext newvc = (new VariantContextBuilder(vc))
                    .attribute("SIMID", "real" + vindex)
                    .make();

            vcw.add(newvc);

            vindex++;
        }

        vcw.close();

        sout.println(dt);
    }
}
