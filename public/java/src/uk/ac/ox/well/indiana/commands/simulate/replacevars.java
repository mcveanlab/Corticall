package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.apache.ivy.util.StringUtils;
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

public class replacevars extends Module {
    @Argument(fullName="ref", shortName="r", doc="Reference sequence")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="refgff", shortName="rg", doc="Reference GFF")
    public GFF3 REFGFF;

    @Argument(fullName="alt", shortName="a", doc="Reference sequence (alt)")
    public IndexedFastaSequenceFile ALT;

    @Argument(fullName="altgff", shortName="ag", doc="Alternate GFF")
    public GFF3 ALTGFF;

    @Argument(fullName="sampleName", shortName="sn", doc="Sample name")
    public String SAMPLE_NAME = "HB3/PG0052-C/ERR019054";

    @Argument(fullName="seed", shortName="s", doc="Seed for RNG", required=false)
    public Long SEED;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 9;

    @Argument(fullName="telomere", shortName="t", doc="Telomere map")
    public File TELOMERE;

    @Output
    public File out;

    @Output(fullName="replacementsOut", shortName="ro", doc="Output file for var replacements")
    public PrintStream rout;

    private List<Interval> loadTelomereMap() {
        List<Interval> tm = new ArrayList<Interval>();

        TableReader tr = new TableReader(TELOMERE);
        for (Map<String, String> te : tr) {
            String chrom = te.get("chrom");
            int start = Integer.valueOf(te.get("co_pos_min"));
            int end = Integer.valueOf(te.get("co_pos_max"));

            Interval interval = new Interval(chrom, start, end);
            tm.add(interval);
        }

        return tm;
    }

    @Override
    public void execute() {
        Map<String, List<GFF3Record>> varRef = new TreeMap<String, List<GFF3Record>>();
        for (GFF3Record refgr : REFGFF) {
            if (refgr.getType().equals("gene")) {
                String upsClass = refgr.getAttribute("class").replaceAll("\\*", "");

                if (upsClass.equals("U") || upsClass.equals("UPSB2") || upsClass.equals("UPSB5")) {
                    upsClass = "UNKNOWN";
                }

                if (!varRef.containsKey(upsClass)) {
                    varRef.put(upsClass, new ArrayList<GFF3Record>());
                }

                varRef.get(upsClass).add(refgr);
            }
        }

        Map<String, List<GFF3Record>> varAlt = new TreeMap<String, List<GFF3Record>>();
        for (GFF3Record altgr : ALTGFF) {
            if (altgr.getType().equals("gene")) {
                String upsClass = altgr.getAttribute("class").replaceAll("\\*", "");

                if (upsClass.equals("ND")) {
                    upsClass = "UNKNOWN";
                }

                if (!varAlt.containsKey(upsClass)) {
                    varAlt.put(upsClass, new ArrayList<GFF3Record>());
                }

                varAlt.get(upsClass).add(altgr);
            }
        }

        Map<GFF3Record, GFF3Record> refToAltMap = new HashMap<GFF3Record, GFF3Record>();

        Map<String, Map<Integer, List<VariantContext>>> variants = new HashMap<String, Map<Integer, List<VariantContext>>>();

        if (SEED == null) { SEED = System.currentTimeMillis(); }
        Random rng = new Random(SEED);

        List<Interval> telomereMap = loadTelomereMap();

        for (String upsClass : varRef.keySet()) {
            List<GFF3Record> refgrs = varRef.get(upsClass);
            List<GFF3Record> altgrs = varAlt.get(upsClass);

            if (refgrs != null && altgrs != null) {
                for (int i = 0; i < (refgrs.size() < altgrs.size() ? refgrs.size() : altgrs.size()); i++) {
                    refToAltMap.put(refgrs.get(i), altgrs.get(i));

                    rout.println(refgrs.get(i).getAttribute("ID") + "\t" + altgrs.get(i).getAttribute("ID"));

                    log.info("  var3D7={} varHB3={}", refgrs.get(i).getAttribute("ID"), altgrs.get(i).getAttribute("ID"));
                }
            }

            if (refgrs != null && refgrs.size() > 1) {
                int index1, index2;
                GFF3Record var1, var2;

                do {
                    index1 = rng.nextInt(refgrs.size());
                    index2 = rng.nextInt(refgrs.size());

                    var1 = refgrs.get(index1);
                    var2 = refgrs.get(index2);
                } while (index1 == index2 || var1.getSeqid().equals(var2.getSeqid()));

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

                    log.info("         nv={}", newVar.toString());
                    log.info("         sw={}", switches.toString());

                    int ti = rng.nextInt(telomereMap.size());
                    Interval interval = telomereMap.get(ti);

                    int pos = rng.nextInt(interval.length()) + interval.getStart();

                    Allele refAllele = Allele.create(new String(REF.getSubsequenceAt(interval.getSequence(), pos, pos).getBases()), true);
                    Allele varAllele = Allele.create(newVar.toString(), false);

                    Genotype insg = (new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(varAllele))).make();
                    GenotypesContext insgc = GenotypesContext.create(insg);

                    VariantContext vci = (new VariantContextBuilder())
                            .chr(interval.getSequence())
                            .start(pos)
                            .stop(pos)
                            .noID()
                            .attribute("VAR", "NAHR_" + var1.getAttribute("ID") + "_" + var2.getAttribute("ID"))
                            .alleles(Arrays.asList(refAllele, varAllele))
                            .genotypes(insgc)
                            .make();

                    if (!variants.containsKey(interval.getSequence())) {
                        variants.put(interval.getSequence(), new TreeMap<Integer, List<VariantContext>>());
                    }

                    if (!variants.get(interval.getSequence()).containsKey(pos)) {
                        variants.get(interval.getSequence()).put(pos, new ArrayList<VariantContext>());
                    }

                    variants.get(interval.getSequence()).get(pos).add(vci);
                }
            }
        }

        for (GFF3Record refgr : refToAltMap.keySet()) {
            GFF3Record altgr = refToAltMap.get(refgr);

            Allele rbef = Allele.create(new String(REF.getSubsequenceAt(refgr.getSeqid(), refgr.getStart() - 1, refgr.getStart() - 1).getBases()));
            Allele rall = Allele.create(new String(REF.getSubsequenceAt(refgr.getSeqid(), refgr.getStart(), refgr.getEnd()).getBases()), true);
            List<Allele> ralls = new ArrayList<Allele>();
            ralls.add(rall);
            ralls.add(rbef);

            Genotype delg = (new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(rbef))).make();
            GenotypesContext delgc = GenotypesContext.create(delg);

            VariantContext vcd = (new VariantContextBuilder())
                    .chr(refgr.getSeqid())
                    .start(refgr.getStart())
                    .computeEndFromAlleles(ralls, refgr.getStart())
                    .noID()
                    .attribute("VAR", refgr.getAttribute("ID"))
                    .alleles(ralls)
                    .genotypes(delgc)
                    .make();

            Allele abef = Allele.create(new String(REF.getSubsequenceAt(refgr.getSeqid(), refgr.getStart() - 1, refgr.getStart() - 1).getBases()), true);
            Allele aall = Allele.create(new String(ALT.getSubsequenceAt(altgr.getSeqid(), altgr.getStart(), altgr.getEnd()).getBases()), false);
            List<Allele> aalls = new ArrayList<Allele>();
            aalls.add(abef);
            aalls.add(aall);

            Genotype insg = (new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(aall))).make();
            GenotypesContext insgc = GenotypesContext.create(insg);

            VariantContext vci = (new VariantContextBuilder())
                    .chr(refgr.getSeqid())
                    .start(refgr.getStart())
                    .stop(refgr.getStart())
                    .noID()
                    .attribute("VAR", altgr.getAttribute("ID"))
                    .alleles(aalls)
                    .genotypes(insgc)
                    .make();

            if (!variants.containsKey(refgr.getSeqid())) {
                variants.put(refgr.getSeqid(), new TreeMap<Integer, List<VariantContext>>());
            }

            if (!variants.get(refgr.getSeqid()).containsKey(refgr.getStart())) {
                variants.get(refgr.getSeqid()).put(refgr.getStart(), new ArrayList<VariantContext>());
            }

            variants.get(refgr.getSeqid()).get(refgr.getStart()).add(vcd);
            variants.get(refgr.getSeqid()).get(refgr.getStart()).add(vci);
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
        header.addMetaDataLine(new VCFInfoHeaderLine("VAR", 1, VCFHeaderLineType.String, "Is a var gene"));

        vcw.writeHeader(header);

        for (SAMSequenceRecord ssr : REF.getSequenceDictionary().getSequences()) {
            String chrom = ssr.getSequenceName();

            if (variants.containsKey(chrom)) {
                for (Integer pos : variants.get(chrom).keySet()) {
                    for (VariantContext vc : variants.get(chrom).get(pos)) {
                        vcw.add(vc);
                    }
                }
            }
        }

        vcw.close();
    }
}
