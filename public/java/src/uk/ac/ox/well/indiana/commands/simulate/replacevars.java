package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

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

    @Output
    public File out;

    @Output(fullName="replacementsOut", shortName="rout", doc="Output file for var replacements")
    public PrintStream rout;

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

        for (String upsClass : varRef.keySet()) {
            List<GFF3Record> refgrs = varRef.get(upsClass);
            List<GFF3Record> altgrs = varAlt.get(upsClass);

            if (refgrs != null && altgrs != null) {
                for (int i = 0; i < (refgrs.size() < altgrs.size() ? refgrs.size() : altgrs.size()); i++) {
                    refToAltMap.put(refgrs.get(i), altgrs.get(i));
                }
            }
        }

        log.info("replacements: {}", refToAltMap.size());

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

        vcw.writeHeader(header);

        Map<String, Map<Integer, List<VariantContext>>> variants = new HashMap<String, Map<Integer, List<VariantContext>>>();

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
                    .attributes(null)
                    .alleles(ralls)
                    .genotypes(delgc)
                    .make();

            //vcw.add(vcd);

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
                    .attributes(null)
                    .alleles(aalls)
                    .genotypes(insgc)
                    .make();

            //vcw.add(vci);

            if (!variants.containsKey(refgr.getSeqid())) {
                variants.put(refgr.getSeqid(), new TreeMap<Integer, List<VariantContext>>());
            }

            if (!variants.get(refgr.getSeqid()).containsKey(refgr.getStart())) {
                variants.get(refgr.getSeqid()).put(refgr.getStart(), new ArrayList<VariantContext>());
            }

            variants.get(refgr.getSeqid()).get(refgr.getStart()).add(vcd);
            variants.get(refgr.getSeqid()).get(refgr.getStart()).add(vci);
        }

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
