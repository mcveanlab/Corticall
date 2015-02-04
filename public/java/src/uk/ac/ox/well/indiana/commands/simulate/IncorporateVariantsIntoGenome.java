package uk.ac.ox.well.indiana.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.apache.ivy.util.StringUtils;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class IncorporateVariantsIntoGenome extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public FastaSequenceFile REF;

    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Output
    public PrintStream out;

    @Output(fullName="vout", shortName="vo", doc="Variants (with flanks)")
    public File vout;

    @Output(fullName="f1out", shortName="f1o", doc="Fastq out 1")
    public PrintStream f1out;

    @Output(fullName="f2out", shortName="f2o", doc="Fastq out 2")
    public PrintStream f2out;

    @Override
    public void execute() {
        Map<String, Map<Integer, Set<VariantContext>>> variants = new HashMap<String, Map<Integer, Set<VariantContext>>>();

        log.info("Loading variants...");
        for (VariantContext vc : VCF) {
            String chr = vc.getChr();
            int start = vc.getStart();

            if (!variants.containsKey(chr)) {
                variants.put(chr, new HashMap<Integer, Set<VariantContext>>());
            }

            if (!variants.get(chr).containsKey(start)) {
                variants.get(chr).put(start - 1, new HashSet<VariantContext>());
            }

            variants.get(chr).get(start - 1).add(vc);
        }

        Map<String, Map<String, String>> attributes = new HashMap<String, Map<String, String>>();

        log.info("Processing bases...");
        ReferenceSequence rseq;
        while ((rseq = REF.nextSequence()) != null) {
            String[] names = rseq.getName().split("\\s+");
            String name = names[0];

            log.info("  {}", name);

            String seq = new String(rseq.getBases());
            List<String> alleles = new ArrayList<String>(seq.length());

            for (int i = 0; i < seq.length(); i++) {
                alleles.add(String.valueOf(seq.charAt(i)));
            }

            for (int i = 0; i < seq.length(); i++) {
                if (variants.containsKey(name) && variants.get(name).containsKey(i)) {
                    for (VariantContext vc : variants.get(name).get(i)) {
                        String alleleToRemove = vc.getReference().getBaseString();
                        String alleleToAdd = !vc.isVariant() ? vc.getReference().getBaseString() : vc.getAlternateAllele(0).getBaseString();

                        for (int j = 0; j < alleleToRemove.length(); j++) {
                            alleles.set(i + j, "");
                        }

                        alleles.set(i, alleleToAdd);
                    }
                }
            }

            for (int i = 0; i < seq.length(); i++) {
                if (variants.containsKey(name) && variants.get(name).containsKey(i)) {
                    for (VariantContext vc : variants.get(name).get(i)) {
                        StringBuilder  leftFlankBuilder = new StringBuilder();
                        StringBuilder rightFlankBuilder = new StringBuilder();

                        int pos = i - 1;
                        while (leftFlankBuilder.length() <= 50 && pos >= 0) {
                            leftFlankBuilder.insert(0, alleles.get(pos));
                            pos--;
                        }

                        pos = i + 1;
                        while (rightFlankBuilder.length() <= 50 && pos < alleles.size()) {
                            rightFlankBuilder.append(alleles.get(pos));
                            pos++;
                        }

                        String id = vc.getAttributeAsString("SIMID", "unknown");

                        if (!attributes.containsKey(id)) {
                            attributes.put(id, new HashMap<String, String>());
                        }

                        String leftFlank = leftFlankBuilder.length() > 50 ? leftFlankBuilder.substring(leftFlankBuilder.length() - 50, leftFlankBuilder.length()) : leftFlankBuilder.toString();
                        String rightFlank = rightFlankBuilder.length() > 50 ? rightFlankBuilder.substring(rightFlankBuilder.length() - 50, rightFlankBuilder.length()) : rightFlankBuilder.toString();

                        attributes.get(id).put("flank_left", leftFlank);
                        attributes.get(id).put("flank_right", rightFlank);
                    }
                }
            }

            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < seq.length(); i++) {
                sb.append(alleles.get(i));
            }

            out.println(">" + name);
            out.println(sb.toString());
        }

        VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder();
        vcwb.setOutputFile(vout);
        vcwb.setReferenceDictionary(VCF.getFileHeader().getSequenceDictionary());
        VariantContextWriter vcw = vcwb.build();

        VCFHeader header = new VCFHeader(VCF.getFileHeader());

        header.addMetaDataLine(new VCFInfoHeaderLine("flank_left", 1, VCFHeaderLineType.String, "The 5' flanking sequence"));
        header.addMetaDataLine(new VCFInfoHeaderLine("flank_right", 1, VCFHeaderLineType.String, "The 3' flanking sequence"));
        header.setSequenceDictionary(VCF.getFileHeader().getSequenceDictionary());

        vcw.writeHeader(header);

        int with = 0, without = 0, variantNum  = 0;
        for (VariantContext vc : VCF) {
            variantNum++;

            String id = vc.getAttributeAsString("SIMID", "unknown");

            String end1id = String.format("@%s variant%d/1", id, variantNum);
            String end2id = String.format("@%s variant%d/2", id, variantNum);

            if (attributes.containsKey(id)) {
                VariantContext newvc = (new VariantContextBuilder(vc))
                        .attribute("flank_left", attributes.get(id).get("flank_left"))
                        .attribute("flank_right", attributes.get(id).get("flank_right"))
                        .make();

                vcw.add(newvc);

                with++;

                f1out.println(end1id);
                f1out.println(attributes.get(id).get("flank_left"));
                f1out.println("+");
                f1out.println(StringUtils.repeat("I", attributes.get(id).get("flank_left").length()));

                f2out.println(end2id);
                f2out.println(attributes.get(id).get("flank_right"));
                f2out.println("+");
                f2out.println(StringUtils.repeat("I", attributes.get(id).get("flank_right").length()));
            } else {
                vcw.add(vc);

                log.info("Missing attributes: {}", vc);

                without++;
            }
        }

        log.info("  with={} without={}", with, without);

        vcw.close();
    }
}
