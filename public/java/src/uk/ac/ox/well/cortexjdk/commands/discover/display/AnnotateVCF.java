package uk.ac.ox.well.cortexjdk.commands.discover.display;

import com.google.common.base.Joiner;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.gff.GFF3;
import uk.ac.ox.well.cortexjdk.utils.io.gff.GFF3Record;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;

import java.io.File;
import java.util.*;

import static htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType.VCF;

public class AnnotateVCF extends Module {
    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VARIANTS;

    @Argument(fullName="core", shortName="c", doc="Core bed")
    public File CORE_BED;

    @Argument(fullName="accessory", shortName="a", doc="Accessory bed")
    public File ACCESSORY_BED;

    @Argument(fullName="genes", shortName="g", doc="Genes")
    public ArrayList<GFF3> GENES;

    @Output
    public File out;

    @Override
    public void execute() {
        IntervalTreeMap<String> itc = new IntervalTreeMap<>();
        IntervalTreeMap<String> ita = new IntervalTreeMap<>();
        IntervalTreeMap<GFF3Record> itg = new IntervalTreeMap<>();

        TableReader trcore = new TableReader(CORE_BED, "chrom", "start", "stop", "label");
        for (Map<String, String> te : trcore) {
            try {
                int x = Integer.valueOf(te.get("start"));
                int y = Integer.valueOf(te.get("stop"));
            } catch (NumberFormatException e) {
                log.info("{}", Joiner.on(" ").withKeyValueSeparator("=").join(te));
            }

            Interval it = new Interval(te.get("chrom"), Integer.valueOf(te.get("start")), Integer.valueOf(te.get("stop")));
            itc.put(it, te.get("label"));
        }

        TableReader tracc = new TableReader(ACCESSORY_BED, "chrom", "start", "stop", "label");
        for (Map<String, String> te : tracc) {
            Interval it = new Interval(te.get("chrom"), Integer.valueOf(te.get("start")), Integer.valueOf(te.get("stop")));
            ita.put(it, te.get("label"));
        }

        for (GFF3 g : GENES) {
            for (GFF3Record gr : g) {
                if (gr.getType().equals("gene")) {
                    Interval it = new Interval(gr.getSeqid(), gr.getStart(), gr.getEnd());
                    itg.put(it, gr);
                }
            }
        }

        VariantContextWriter vcw = new VariantContextWriterBuilder()
            .setOutputFile(out)
            .setOutputFileType(VCF)
            .setOption(Options.DO_NOT_WRITE_GENOTYPES)
            .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
            .unsetOption(Options.INDEX_ON_THE_FLY)
            .build();

        vcw.writeHeader(VARIANTS.getFileHeader());

        for (VariantContext vc : VARIANTS) {
            Interval it = new Interval(vc.getContig(), vc.getStart(), vc.getEnd());

            String label = "unknown";
            if (itc.containsOverlapping(it)) { label = "core"; }
            if (ita.containsOverlapping(it)) { label = "accessory"; }

            Set<String> genes = new TreeSet<>();
            //Set<String> desc = new TreeSet<>();
            for (GFF3Record gr : itg.getOverlapping(it)) {
                genes.add(gr.getAttribute("ID"));
                //desc.add(gr.getAttribute("description"));
            }

            VariantContext newvc = new VariantContextBuilder(vc)
                    .attribute("REGION", label)
                    .attribute("GENES", Joiner.on(",").join(genes))
                    //.attribute("DESC", Joiner.on(",").join(desc))
                    .make();

            vcw.add(newvc);
        }

        vcw.close();
    }
}
