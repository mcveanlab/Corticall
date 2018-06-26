package uk.ac.ox.well.cortexjdk.commands.discover.call;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
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
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;

import java.io.File;
import java.util.*;

import static htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType.VCF;

public class AnnotateCalls extends Module {
    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VARIANTS;

    @Argument(fullName="core", shortName="c", doc="Core bed")
    public File CORE_BED;

    @Argument(fullName="accessory", shortName="a", doc="Accessory bed")
    public File ACCESSORY_BED;

    @Argument(fullName="genes", shortName="g", doc="Genes")
    public ArrayList<GFF3> GENES;

    @Argument(fullName="partitions", shortName="p", doc="Partitions")
    public FastaSequenceFile PARTITIONS;

    @Argument(fullName="rois", shortName="r", doc="ROIs")
    public CortexGraph ROIS;

    @Output
    public File out;

    @Override
    public void execute() {
        Set<CanonicalKmer> cks = new HashSet<>();
        for (CortexRecord cr : ROIS) {
            cks.add(cr.getCanonicalKmer());
        }

        Map<String, String> refs = new HashMap<>();
        ReferenceSequence rseq;
        while ((rseq = PARTITIONS.nextSequence()) != null) {
            refs.put(rseq.getName().split(" ")[0], rseq.getBaseString());
        }

        IntervalTreeMap<String> itc = new IntervalTreeMap<>();
        IntervalTreeMap<String> ita = new IntervalTreeMap<>();
        IntervalTreeMap<GFF3Record> itg = new IntervalTreeMap<>();

        TableReader trcore = new TableReader(CORE_BED, "chrom", "start", "stop", "label");
        for (Map<String, String> te : trcore) {
//            try {
//                int x = Integer.valueOf(te.get("start"));
//                int y = Integer.valueOf(te.get("stop"));
//            } catch (NumberFormatException e) {
//                log.info("{}", Joiner.on(" ").withKeyValueSeparator("=").join(te));
//            }

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
                if (gr.getType().contains("gene")) {
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

            Set<String> closest = new TreeSet<>();
            Interval itn = new Interval(vc.getContig(), vc.getStart() - 100000, vc.getEnd() + 100000);
            List<GFF3Record> grs = new ArrayList<>(itg.getOverlapping(itn));
            if (grs.size() > 0) {
                grs.sort((g1, g2) -> {
                    if (g1.getStart() == g2.getStart()) { return 0; }
                    return Math.abs(g1.getStart() - vc.getStart()) < Math.abs(g2.getStart() - vc.getStart()) ? -1 : 1;
                });

                closest.add(grs.get(0).getAttribute("ID"));
            }

            String pname = vc.getAttributeAsString("PARTITION_NAME", "");
            int plength = 0;
            int numNovels = 0;
            if (refs.containsKey(pname)) {
                plength = refs.get(pname).length();

                for (int i = 0; i <= refs.get(pname).length() - ROIS.getKmerSize(); i++) {
                    CanonicalKmer ck = new CanonicalKmer(refs.get(pname).substring(i, i + ROIS.getKmerSize()));

                    if (cks.contains(ck)) {
                        numNovels++;
                    }
                }
            }

            VariantContext newvc = new VariantContextBuilder(vc)
                    .attribute("REGION", label)
                    .attribute("GENE", Joiner.on(",").join(genes))
                    .attribute("CLOSEST_GENE", Joiner.on(",").join(closest))
                    .attribute("PARTITION_LENGTH", plength)
                    .attribute("PARTITION_NOVELS", numNovels)
                    //.attribute("DESC", Joiner.on(",").join(desc))
                    .make();

            vcw.add(newvc);
        }

        vcw.close();
    }
}
