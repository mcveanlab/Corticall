package uk.ac.ox.well.indiana.attic.analyses.nahr;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.PrintStream;
import java.util.LinkedHashMap;
import java.util.Map;

public class PaintTransmittedHaplotypes extends Module {
    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Argument(fullName="parent1", shortName="p1", doc="Parent 1 sample name")
    public String PARENT1;

    @Argument(fullName="parent2", shortName="p2", doc="Parent 2 sample name")
    public String PARENT2;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableWriter tw = new TableWriter(out);

        log.info("Processing variants...");
        int numVariants = 0;
        for (VariantContext vc : VCF) {
            if (numVariants % 10000 == 0) {
                log.info("  {}", numVariants);
            }
            numVariants++;

            if (!vc.isFiltered()) {
                Genotype g1 = vc.getGenotype(PARENT1);
                Genotype g2 = vc.getGenotype(PARENT2);

                Map<String, String> te = new LinkedHashMap<String, String>();

                te.put("CHROM", vc.getChr());
                te.put("POS", String.valueOf(vc.getStart()));

                for (String sampleName : vc.getSampleNames()) {
                    if (sampleName.equals(PARENT1)) {
                        te.put(PARENT1, "1");
                    } else if (sampleName.equals(PARENT2)) {
                        te.put(PARENT2, "2");
                    } else {
                        Genotype gp = vc.getGenotype(sampleName);

                        int parentIndex = 4;

                        if ( g1.isNoCall() || g2.isNoCall() || gp.isNoCall() ) {
                            parentIndex = 0;
                        } else if ( ( (gp.getPloidy() == 1) || (g1.isPhased() && g2.isPhased() && gp.isPhased()) ) ) {
                            if (gp.sameGenotype(g1, false) && !gp.sameGenotype(g2, false)) {
                                parentIndex = 1;
                            } else if (gp.sameGenotype(g2, false) && !gp.sameGenotype(g1, false)) {
                                parentIndex = 2;
                            }
                        }

                        te.put(sampleName, String.valueOf(parentIndex));
                    }
                }

                tw.addEntry(te);
            }
        }
    }
}
