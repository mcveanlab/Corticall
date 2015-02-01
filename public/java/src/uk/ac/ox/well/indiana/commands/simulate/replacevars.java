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

public class ReplaceVars extends Module {
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

    @Output(fullName="replacementsOut", shortName="ro", doc="Output file for var replacements")
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

        Map<String, Map<Integer, List<VariantContext>>> variants = new HashMap<String, Map<Integer, List<VariantContext>>>();

        for (String upsClass : varRef.keySet()) {
            List<GFF3Record> refgrs = varRef.get(upsClass);
            List<GFF3Record> altgrs = varAlt.get(upsClass);

            if (refgrs != null && altgrs != null) {
                for (int i = 0; i < (refgrs.size() < altgrs.size() ? refgrs.size() : altgrs.size()); i++) {
                    refToAltMap.put(refgrs.get(i), altgrs.get(i));

                    rout.println(refgrs.get(i).getAttribute("ID") + "\t" + refgrs.get(i).getAttribute("class") + "\t" + altgrs.get(i).getAttribute("ID") + "\t" + altgrs.get(i).getAttribute("class"));

                    log.info("  var3D7={} varHB3={}", refgrs.get(i).getAttribute("ID"), altgrs.get(i).getAttribute("ID"));
                }
            }
        }

        for (GFF3Record refgr : refToAltMap.keySet()) {
            GFF3Record altgr = refToAltMap.get(refgr);

            /*
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
            */

            Allele oldVar = Allele.create(new String(REF.getSubsequenceAt(refgr.getSeqid(), refgr.getStart(), refgr.getEnd()).getBases()), true);
            Allele newVar = Allele.create(new String(ALT.getSubsequenceAt(altgr.getSeqid(), altgr.getStart(), altgr.getEnd()).getBases()), false);

            Genotype varg = (new GenotypeBuilder(SAMPLE_NAME, Arrays.asList(newVar))).make();

            VariantContext vcv = (new VariantContextBuilder())
                    .chr(refgr.getSeqid())
                    .start(refgr.getStart())
                    .computeEndFromAlleles(Arrays.asList(oldVar, newVar), refgr.getStart())
                    .noID()
                    .attribute("VAR_OLD", refgr.getAttribute("ID"))
                    .attribute("VAR_OLD_CLASS", refgr.getAttribute("class"))
                    .attribute("VAR_OLD_POS", refgr.getAttribute("position"))
                    .attribute("VAR_NEW", altgr.getAttribute("ID"))
                    .attribute("VAR_NEW_CLASS", altgr.getAttribute("class"))
                    .attribute("VAR_NEW_POS", altgr.getAttribute("position"))
                    .alleles(Arrays.asList(oldVar, newVar))
                    .genotypes(varg)
                    .make();

            if (!variants.containsKey(refgr.getSeqid())) {
                variants.put(refgr.getSeqid(), new TreeMap<Integer, List<VariantContext>>());
            }

            if (!variants.get(refgr.getSeqid()).containsKey(refgr.getStart())) {
                variants.get(refgr.getSeqid()).put(refgr.getStart(), new ArrayList<VariantContext>());
            }

            variants.get(refgr.getSeqid()).get(refgr.getStart()).add(vcv);
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
        header.addMetaDataLine(new VCFInfoHeaderLine("VAR_OLD", 1, VCFHeaderLineType.String, "The var gene that was removed"));
        header.addMetaDataLine(new VCFInfoHeaderLine("VAR_OLD_CLASS", 1, VCFHeaderLineType.String, "The old var UPS class"));
        header.addMetaDataLine(new VCFInfoHeaderLine("VAR_OLD_POS", 1, VCFHeaderLineType.String, "The old var position"));
        header.addMetaDataLine(new VCFInfoHeaderLine("VAR_NEW", 1, VCFHeaderLineType.String, "The var gene replacement"));
        header.addMetaDataLine(new VCFInfoHeaderLine("VAR_NEW_CLASS", 1, VCFHeaderLineType.String, "The new var UPS class"));
        header.addMetaDataLine(new VCFInfoHeaderLine("VAR_NEW_POS", 1, VCFHeaderLineType.String, "The new var position"));

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
