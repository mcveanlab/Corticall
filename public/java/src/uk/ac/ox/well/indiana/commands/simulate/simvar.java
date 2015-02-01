package uk.ac.ox.well.indiana.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang.StringUtils;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class SimVar extends Module {
    @Argument(fullName="ref", shortName="r", doc="Reference sequence")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="refgff", shortName="rg", doc="Reference GFF")
    public GFF3 REFGFF;

    @Argument(fullName="seed", shortName="s", doc="Seed for RNG")
    public Long SEED = System.currentTimeMillis();

    @Argument(fullName="sampleName", shortName="sn", doc="Sample name")
    public String SAMPLE_NAME = "HB3/PG0052-C/ERR019054";

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 9;

    @Argument(fullName="vcf", shortName="v", doc="Child VCF")
    public VCFFileReader VCF;

    @Output
    public File out;

    @Output(fullName="sout", shortName="so", doc="Sequence recomb out")
    public PrintStream sout;

    @Override
    public void execute() {
        Random rng = new Random(SEED);

        Map<String, GFF3Record> refVarIDs = new HashMap<String, GFF3Record>();
        for (GFF3Record gr : REFGFF) {
            if (gr.getType().equals("gene")) {
                refVarIDs.put(gr.getAttribute("ID"), gr);
            }
        }

        for (VariantContext vc : VCF) {
            if (vc.hasAttribute("VAR_OLD")) {
                String oldVarId = vc.getAttributeAsString("VAR_OLD", "");

                refVarIDs.remove(oldVarId);
            }
        }

        Map<String, Map<String, List<GFF3Record>>> similarVars = new HashMap<String, Map<String, List<GFF3Record>>>();

        for (String id : refVarIDs.keySet()) {
            GFF3Record gr = refVarIDs.get(id);
            String upsClass = gr.getAttribute("class");
            String pos = gr.getAttribute("position");

            if (!similarVars.containsKey(upsClass)) { similarVars.put(upsClass, new HashMap<String, List<GFF3Record>>()); }
            if (!similarVars.get(upsClass).containsKey(pos)) { similarVars.get(upsClass).put(pos, new ArrayList<GFF3Record>()); }

            if (pos.contains("telomere")) {
                similarVars.get(upsClass).get(pos).add(gr);
            }
        }

        sout.println(Joiner.on("\t").join(
                "id1",
                "id2",
                "class1",
                "class2",
                "position1",
                "position2",
                "seq1",
                "seq2",
                "sites1",
                "sites2",
                "recombs"
        ));

        Map<String, Map<Integer, List<VariantContext>>> variants = new HashMap<String, Map<Integer, List<VariantContext>>>();

        for (String upsClass : similarVars.keySet()) {
            for (String pos : similarVars.get(upsClass).keySet()) {
                if (similarVars.get(upsClass).get(pos).size() >= 2) {
                    List<GFF3Record> varList = similarVars.get(upsClass).get(pos);

                    GFF3Record var1, var2;

                    do {
                        int index1 = rng.nextInt(varList.size());
                        int index2 = rng.nextInt(varList.size());

                        var1 = varList.get(index1);
                        var2 = varList.get(index2);
                    } while (var1.getSeqid().equals(var2.getSeqid()));

                    log.info("  recombine {}", var1);
                    log.info("            {}", var2);

                    String seq1 = new String(REF.getSubsequenceAt(var1.getSeqid(), var1.getStart(), var1.getEnd()).getBases());
                    String seq2 = new String(REF.getSubsequenceAt(var2.getSeqid(), var2.getStart(), var2.getEnd()).getBases());

                    if (var1.getStrand() != var2.getStrand()) {
                        seq2 = SequenceUtils.reverseComplement(seq1);
                    }

                    log.info("  recombine {} ({}) {} ({})", var1.getAttribute("ID"), var1.getStrand(), var2.getAttribute("ID"), var2.getStrand());
                    log.info("       var1={}", seq1);
                    log.info("       var2={}", seq2);

                    Map<String, Integer> kmers1 = new HashMap<String, Integer>();
                    for (int i = 0; i <= seq1.length() - KMER_SIZE; i++) {
                        String kmer = seq1.substring(i, i + KMER_SIZE);

                        if (!kmers1.containsKey(kmer)) {
                            kmers1.put(kmer, i);
                        }
                    }

                    Map<String, Integer> kmers2 = new HashMap<String, Integer>();
                    for (int i = 0; i <= seq2.length() - KMER_SIZE; i++) {
                        String kmer = seq2.substring(i, i + KMER_SIZE);

                        if (!kmers2.containsKey(kmer)) {
                            kmers2.put(kmer, i);
                        }
                    }

                    List<Integer> possibleRecombSites1 = new ArrayList<Integer>();
                    List<Integer> possibleRecombSites2 = new ArrayList<Integer>();

                    StringBuilder sb = new StringBuilder();
                    for (int i = 0; i <= (seq1.length() < seq2.length() ? seq1.length() : seq2.length()) - KMER_SIZE; i++) {
                        String kmer1 = seq1.substring(i, i + KMER_SIZE);

                        sb.append(kmers1.containsKey(kmer1) && kmers2.containsKey(kmer1) ? "1" : "0");

                        if (kmers1.containsKey(kmer1) && kmers2.containsKey(kmer1)) {
                            int pos1 = i;
                            int pos2 = kmers2.get(kmer1);

                            if (Math.abs(pos1 - pos2) < 100 && Math.abs(pos1 - pos2) > 20) {
                                possibleRecombSites1.add(pos1);
                                possibleRecombSites2.add(pos2);
                            }
                        }
                    }

                    log.info("         sb={}", sb.toString());

                    int numRecombs = rng.nextInt(3) + 2;

                    if (possibleRecombSites1.size() > numRecombs) {
                        Set<Integer> recombs = new TreeSet<Integer>();
                        for (int nr = 0; nr < numRecombs; nr++) {
                            recombs.add(rng.nextInt(possibleRecombSites1.size()));
                        }

                        boolean copyFrom1 = true;
                        int lastPos = 0;

                        StringBuilder newVar = new StringBuilder();
                        StringBuilder switches = new StringBuilder();
                        for (Integer index : recombs) {
                            if (copyFrom1) {
                                int currentPos = possibleRecombSites1.get(index);
                                newVar.append(seq1.substring(lastPos, currentPos));
                                switches.append(StringUtils.repeat(" ", currentPos - lastPos - 1)).append("v");

                                lastPos = possibleRecombSites2.get(index);
                            } else {
                                int currentPos = possibleRecombSites2.get(index);
                                newVar.append(seq2.substring(lastPos, currentPos));
                                switches.append(StringUtils.repeat(" ", currentPos - lastPos - 1)).append("^");

                                lastPos = possibleRecombSites1.get(index);
                            }

                            copyFrom1 = !copyFrom1;
                        }

                        if (copyFrom1) {
                            newVar.append(seq1.substring(lastPos, seq1.length()));
                        } else {
                            newVar.append(seq2.substring(lastPos, seq2.length()));
                        }

                        sout.println(Joiner.on("\t").join(
                                var1.getAttribute("ID"),
                                var2.getAttribute("ID"),
                                var1.getAttribute("class"),
                                var2.getAttribute("class"),
                                var1.getAttribute("position"),
                                var2.getAttribute("position"),
                                seq1,
                                seq2,
                                Joiner.on(",").join(possibleRecombSites1),
                                Joiner.on(",").join(possibleRecombSites2),
                                Joiner.on(",").join(recombs)
                        ));

                        log.info("         nv={}", newVar.toString());
                        log.info("         sw={}", switches.toString());

                        Allele oldVar1 = Allele.create(seq1, true);
                        Allele newVar1 = Allele.create(newVar.toString(), false);

                        Genotype newVar1G = (new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(newVar1))).make();

                        VariantContext newVar1VC = (new VariantContextBuilder())
                                .chr(var1.getSeqid())
                                .start(var1.getStart())
                                .computeEndFromAlleles(Arrays.asList(oldVar1, newVar1), var1.getStart())
                                .noID()
                                .attribute("NAHR", var1.getAttribute("ID") + "_" + var2.getAttribute("ID"))
                                .alleles(Arrays.asList(oldVar1, newVar1))
                                .genotypes(newVar1G)
                                .make();

                        String refBase = new String(REF.getSubsequenceAt(var2.getSeqid(), var2.getStart() - 1, var2.getStart() - 1).getBases());

                        Allele oldVar2 = Allele.create(refBase + seq2, true);
                        Allele newVar2 = Allele.create(refBase, false);

                        Genotype newVar2G = (new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(newVar2))).make();

                        VariantContext newVar2VC = (new VariantContextBuilder())
                                .chr(var2.getSeqid())
                                .start(var2.getStart())
                                .computeEndFromAlleles(Arrays.asList(oldVar2, newVar2), var2.getStart())
                                .noID()
                                .attribute("NAHR", var1.getAttribute("ID") + "_" + var2.getAttribute("ID"))
                                .alleles(Arrays.asList(oldVar2, newVar2))
                                .genotypes(newVar2G)
                                .make();

                        if (!variants.containsKey(newVar1VC.getChr())) {
                            variants.put(newVar1VC.getChr(), new TreeMap<Integer, List<VariantContext>>());
                        }

                        if (!variants.get(newVar1VC.getChr()).containsKey(newVar1VC.getStart())) {
                            variants.get(newVar1VC.getChr()).put(newVar1VC.getStart(), new ArrayList<VariantContext>());
                        }

                        variants.get(newVar1VC.getChr()).get(newVar1VC.getStart()).add(newVar1VC);

                        if (!variants.containsKey(newVar2VC.getChr())) {
                            variants.put(newVar2VC.getChr(), new TreeMap<Integer, List<VariantContext>>());
                        }

                        if (!variants.get(newVar2VC.getChr()).containsKey(newVar2VC.getStart())) {
                            variants.get(newVar2VC.getChr()).put(newVar2VC.getStart(), new ArrayList<VariantContext>());
                        }

                        variants.get(newVar2VC.getChr()).get(newVar2VC.getStart()).add(newVar2VC);
                    }
                }
            }
        }

        VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder();
        vcwb.setOutputFile(out);
        vcwb.setReferenceDictionary(REF.getSequenceDictionary());
        VariantContextWriter vcw = vcwb.build();

        Set<String> sampleNames = new HashSet<String>();
        sampleNames.add(SAMPLE_NAME);

        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();

        VCFHeader header = new VCFHeader(headerLines, sampleNames);
        header.setSequenceDictionary(REF.getSequenceDictionary());
        header.addMetaDataLine(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
        header.addMetaDataLine(new VCFInfoHeaderLine("NAHR", 1, VCFHeaderLineType.String, "Is a var recombination"));

        vcw.writeHeader(header);

        for (SAMSequenceRecord ssr : REF.getSequenceDictionary().getSequences()) {
            String chr = ssr.getSequenceName();
            if (variants.containsKey(chr)) {
                for (int pos : variants.get(chr).keySet()) {
                    for (VariantContext vc : variants.get(chr).get(pos)) {
                        vcw.add(vc);
                    }
                }
            }
        }

        vcw.close();
    }
}
