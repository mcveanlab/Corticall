package uk.ac.ox.well.indiana.commands.evaluate;

import com.google.common.base.Joiner;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.SortingVariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.io.File;
import java.io.PrintStream;
import java.util.Map;

public class SelectDeNovosFromVCF extends Module {
    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="mother", shortName="m", doc="Mother")
    public String MOTHER;

    @Argument(fullName="father", shortName="f", doc="Father")
    public String FATHER;

    @Argument(fullName="regionLabels", shortName="l", doc="Region labels")
    public File REGIONS;

    @Output
    public File out;

    private IntervalTreeMap<String> loadRegions() {
        TableReader tr = new TableReader(REGIONS, "chr", "start", "stop", "region");

        IntervalTreeMap<String> itm = new IntervalTreeMap<String>();

        for (Map<String, String> te : tr) {
            String chr = te.get("chr");
            int start = Integer.valueOf(te.get("start"));
            int stop = Integer.valueOf(te.get("stop"));
            String region = te.get("region");

            Interval it = new Interval(chr, start, stop);

            itm.put(it, region);
        }

        return itm;
    }

    @Override
    public void execute() {
        VariantContextWriter vcw = new VariantContextWriterBuilder()
                .setOutputFile(out)
                .clearIndexCreator()
                .setReferenceDictionary(VCF.getFileHeader().getSequenceDictionary())
                .build();

        VCFHeader vh = VCF.getFileHeader();
        vh.addMetaDataLine(new VCFInfoHeaderLine("region", 1, VCFHeaderLineType.String, "region type in which the variant occurs"));

        vcw.writeHeader(VCF.getFileHeader());

        IntervalTreeMap<String> itm = loadRegions();

        for (VariantContext vc : VCF) {
            if ( vc.getGenotype(CHILD).isCalled() && vc.getGenotype(MOTHER).isCalled() && vc.getGenotype(FATHER).isCalled() &&
                !vc.getGenotype(CHILD).getType().equals(vc.getGenotype(MOTHER).getType()) &&
                !vc.getGenotype(CHILD).getType().equals(vc.getGenotype(FATHER).getType()) &&
                 vc.getGenotype(CHILD).getGQ() > 90 && vc.getGenotype(MOTHER).getGQ() > 90 && vc.getGenotype(FATHER).getGQ() > 90 &&
                 vc.getGenotype(CHILD).getDP() > 30 && vc.getGenotype(MOTHER).getDP() > 30 && vc.getGenotype(FATHER).getDP() > 30
                ) {
                Interval it = new Interval(vc.getChr(), vc.getStart(), vc.getEnd());

                int[] adc = vc.getGenotype(CHILD).getAD();
                int[] adm = vc.getGenotype(MOTHER).getAD();
                int[] adf = vc.getGenotype(FATHER).getAD();

                boolean bigdiff = false;

                if (adc.length == 1) {
                    bigdiff = true;
                } else if (adc.length > 1) {
                    if (adc[0] < adc[1]/10 || adc[1] < adc[0]/10) {
                        bigdiff = true;
                    }
                }

                if (bigdiff) {
                    String region = "unknown";

                    if (itm.containsOverlapping(it)) {
                        region = Joiner.on(",").join(itm.getOverlapping(it));
                    } else if (itm.containsContained(it)) {
                        region = Joiner.on(",").join(itm.getContained(it));
                    }

                    VariantContext newvc = new VariantContextBuilder(vc).attribute("region", region).make();

                    vcw.add(newvc);
                }
            }
        }

        vcw.close();
    }
}
