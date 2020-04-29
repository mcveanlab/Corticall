package uk.ac.ox.well.cortexjdk.commands.discover.call;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
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
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.File;
import java.util.*;

import static htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType.VCF;

public class AnnotateCalls extends Module {
    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VARIANTS;

    @Argument(fullName="accessory", shortName="a", doc="Accessory bed")
    public File ACCESSORY_BED;

    @Argument(fullName="repeatmasks", shortName="rm", doc="Repeat masks")
    public ArrayList<GFF3> REPEAT_MASKS;

    @Argument(fullName="genes", shortName="g", doc="Genes")
    public ArrayList<GFF3> GENES;

    @Argument(fullName="partitions", shortName="p", doc="Partitions")
    public FastaSequenceFile PARTITIONS;

    @Argument(fullName="rois", shortName="r", doc="ROIs")
    public CortexGraph ROIS;

    @Argument(fullName="reference", shortName="R", doc="Reference")
    public ArrayList<IndexedFastaSequenceFile> REFERENCE;

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
        IntervalTreeMap<GFF3Record> ite = new IntervalTreeMap<>();

        /*
        TableReader trcore = new TableReader(CORE_BED, "chrom", "start", "stop", "label");
        for (Map<String, String> te : trcore) {
            Interval it = new Interval(te.get("chrom"), Integer.valueOf(te.get("start")), Integer.valueOf(te.get("stop")));
            itc.put(it, te.get("label"));
        }
        */

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

                if (gr.getType().contains("exon")) {
                    Interval it = new Interval(gr.getSeqid(), gr.getStart(), gr.getEnd());
                    ite.put(it, gr);
                }
            }
        }

        IntervalTreeMap<String> itrm = new IntervalTreeMap<>();
        for (GFF3 rm : REPEAT_MASKS) {
            for (GFF3Record r : rm) {
                Interval it = new Interval(r.getSeqid(), r.getStart(), r.getEnd());
                itrm.put(it, r.getAttribute("ID"));
            }
        }

        log.info("Loaded {} repeats", itrm.size());

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

            String label = "core";
            //if (itc.containsOverlapping(it)) { label = "core"; }
            if (ita.containsOverlapping(it)) { label = "accessory"; }

            Set<String> genes = new TreeSet<>();
            for (GFF3Record gr : itg.getOverlapping(it)) {
                genes.add(gr.getAttribute("ID"));

                for (GFF3 g : GENES) {
                    List<GFF3Record> gs = new ArrayList<>(GFF3.getType("exon", g.getContained(gr)));
                    if (!gs.isEmpty()) {
                        gs.sort((g1, g2) -> {
                            if (g1.getStart() == g2.getStart()) { return 0; }
                            return g1.getStart() < g2.getStart() ? -1 : 1;
                        });

                        if (!vc.isSymbolic()) {
                            for (IndexedFastaSequenceFile ref : REFERENCE) {
                                for (SAMSequenceRecord ssra : ref.getSequenceDictionary().getSequences()) {
                                    if (ssra.getSequenceName().equals(gr.getSeqid())) {
                                        List<List<String>> chr = new ArrayList<>();
                                        String seq = ref.getSequence(gr.getSeqid()).getBaseString();
                                        for (int i = 0; i < seq.length(); i++) {
                                            List<String> ns = new ArrayList<>();
                                            ns.add(String.valueOf(seq.charAt(i)).toUpperCase());
                                            chr.add(ns);
                                        }

                                        StringBuilder ob = new StringBuilder();
                                        for (GFF3Record exongr : gs) {
                                            for (int i = exongr.getStart() - 1; i < exongr.getEnd(); i++) {
                                                ob.append(chr.get(i).get(0));
                                            }
                                        }

                                        if (gr.getStrand() == GFF3Record.Strand.NEGATIVE) {
                                            StringBuilder rb = new StringBuilder(SequenceUtils.reverseComplement(ob.toString()));
                                            ob = rb;
                                        }

                                        for (int i = vc.getStart() - 1; i < vc.getEnd(); i++) {
                                            chr.set(i, Arrays.asList(""));
                                        }
                                        chr.set(vc.getStart() - 1, Arrays.asList(vc.getAlternateAllele(0).getBaseString().toLowerCase()));

                                        StringBuilder sb = new StringBuilder();
                                        for (GFF3Record exongr : gs) {
                                            for (int i = exongr.getStart() - 1; i < exongr.getEnd(); i++) {
                                                sb.append(chr.get(i).get(0));
                                            }
                                        }

                                        if (gr.getStrand() == GFF3Record.Strand.NEGATIVE) {
                                            StringBuilder rb = new StringBuilder(SequenceUtils.reverseComplement(sb.toString()));
                                            sb = rb;
                                        }

                                        log.info(ob.toString());
                                        log.info(sb.toString());

                                        int ca = 0;
                                        int cc = 0;
                                        for (int i = 0; i < ob.length(); i++, cc++) {
                                            if (cc > 2) { ca++; cc = 0; }
                                            if (sb.charAt(i) != ob.charAt(i)) {
                                                log.info("{} {} {} {} {}", i, ca, cc, sb.charAt(i), ob.charAt(i));
                                            }
                                        }

                                        log.info("");
                                    }
                                }
                            }
                        }
                    }
                }
            }

            Set<String> closestGene = new TreeSet<>();
            Interval itn = new Interval(vc.getContig(), vc.getStart() - 100000, vc.getEnd() + 100000);
            List<GFF3Record> grs = new ArrayList<>(itg.getOverlapping(itn));
            if (grs.size() > 0) {
                grs.sort((g1, g2) -> {
                    if (g1.getStart() == g2.getStart()) { return 0; }
                    return Math.abs(g1.getStart() - vc.getStart()) < Math.abs(g2.getStart() - vc.getStart()) ? -1 : 1;
                });

                closestGene.add(grs.get(0).getAttribute("ID"));
            }

            String repeat = "NA";
            List<String> repeats = new ArrayList<>(itrm.getOverlapping(it));
            if (repeats.size() > 0) {
                repeat = repeats.get(0);
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
                    .attribute("GENIC", genes.size() > 0 ? "true" : "false")
                    .attribute("GENE", Joiner.on(",").join(genes))
                    .attribute("CLOSEST_GENE", Joiner.on(",").join(closestGene))
                    .attribute("REPEAT", repeat)
                    .attribute("PARTITION_LENGTH", plength)
                    .attribute("PARTITION_NOVELS", numNovels)
                    .make();

            vcw.add(newvc);
        }

        vcw.close();
    }
}
