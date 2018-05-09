package uk.ac.ox.well.cortexjdk.commands.simulate;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

import static htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType.VCF;

public class SimToVCF extends Module {
    @Argument(fullName="sim", shortName="s", doc="Simulation")
    public File SIM;

    @Argument(fullName="background", shortName="b", doc="Background", required=false)
    public TreeMap<Integer, IndexedReference> BACKGROUNDS;

    @Argument(fullName="variants", shortName="v", doc="Variants")
    public VCFFileReader VARIANTS;

    @Output
    public File out;

    @Override
    public void execute() {
        VariantContextWriter vcw = new VariantContextWriterBuilder()
                .setOutputFile(out)
                .setOutputFileType(VCF)
                .setOption(Options.DO_NOT_WRITE_GENOTYPES)
                .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .unsetOption(Options.INDEX_ON_THE_FLY)
                .build();

        vcw.writeHeader(VARIANTS.getFileHeader());

        Map<String, Integer> sid = new HashMap<>();
        for (int i = 0; i < VARIANTS.getFileHeader().getSequenceDictionary().size(); i++) {
            sid.put(VARIANTS.getFileHeader().getSequenceDictionary().getSequence(i).getSequenceName(), i);
        }

        Set<VariantContext> svcs = new TreeSet<>((v1, v2) -> {
            int sid0 = sid.get(v1.getContig());
            int sid1 = sid.get(v2.getContig());

            if (sid0 != sid1) { return sid0 < sid1 ? -1 : 1; }
            if (v1.getStart() != v2.getStart()) { return v1.getStart() < v2.getStart() ? -1 : 1; }

            return 0;
        });

        TableReader tr = new TableReader(SIM);

        for (Map<String, String> te : tr) {
            if (!te.get("type").equals("RECOMB")) {
                log.info("{}", te);

                String sleft = te.get("sleft");
                String sright = te.get("sright");
                String oldAllele = te.get("old").replaceAll("\\.", "");
                String newAllele = te.get("new").replaceAll("\\.", "");

                String oldHap = sleft + oldAllele + sright;
                String newHap = sleft + newAllele + sright;

                IndexedReference ref = BACKGROUNDS.get(Integer.valueOf(te.get("parent")));
                List<SAMRecord> srs = sort(ref.align(oldHap.toUpperCase()));
                SAMRecord srBest = srs.size() > 0 ? srs.iterator().next() : null;

                if (srBest != null) {
                    log.info("{}", sleft);
                    log.info("{}{}", StringUtil.repeatCharNTimes(' ', sleft.length()), oldAllele);
                    log.info("{}{}", StringUtil.repeatCharNTimes(' ', sleft.length() + oldAllele.length()), sright);

                    log.info("{}", sleft);
                    log.info("{}{}", StringUtil.repeatCharNTimes(' ', sleft.length()), newAllele);
                    log.info("{}{}", StringUtil.repeatCharNTimes(' ', sleft.length() + newAllele.length()), sright);

                    int pos = te.get("type").equalsIgnoreCase("SNV") ? srBest.getAlignmentStart() + sleft.length() + 1 : srBest.getAlignmentStart() + sleft.length();

                    String expRefBase = te.get("type").equals("SNV") ? "" : sleft.substring(sleft.length() - 1, sleft.length()).toUpperCase();
                    String actRefBase = (te.get("type").equals("SNV")) ? expRefBase : ref.getReferenceSequence().getSubsequenceAt(srBest.getContig(), pos, pos).getBaseString();

                    log.info("{}{}", StringUtil.repeatCharNTimes(' ', sleft.length() - 1), expRefBase);
                    log.info("{}{}", StringUtil.repeatCharNTimes(' ', sleft.length() - 1), actRefBase);

                    if (!oldAllele.equals(newAllele)) {
                        List<Allele> alleles = Arrays.asList(Allele.create(actRefBase + oldAllele, true), Allele.create(actRefBase + newAllele));

                        VariantContext vc = new VariantContextBuilder()
                                .chr(srBest.getContig())
                                .alleles(alleles)
                                .start(pos)
                                .computeEndFromAlleles(alleles, pos)
                                .attribute("SIM_TYPE", te.get("type"))
                                .attribute("INDEX", te.get("index"))
                                .attribute("OLD_HAP", oldHap.toUpperCase())
                                .attribute("NEW_HAP", newHap.toUpperCase())
                                .attribute("SLEFT", sleft.toUpperCase())
                                .attribute("SOLD", oldAllele)
                                .attribute("SNEW", newAllele)
                                .attribute("SRIGHT", sright.toUpperCase())
                                .noGenotypes()
                                .make();

                        log.info("{}", vc);
                        log.info("");

                        svcs.add(vc);
                    }
                } else {
                    log.info("skipped");
                }
            }
        }

        for (VariantContext vc : svcs) {
            vcw.add(vc);
        }

        vcw.close();
    }

    private List<SAMRecord> sort(List<SAMRecord> srs) {
        srs.sort((o1, o2) -> {
            if (o1.getMappingQuality() != o2.getMappingQuality()) {
                return o1.getMappingQuality() > o2.getMappingQuality() ? -1 : 1;
            }

            if (!o1.getIntegerAttribute("NM").equals(o2.getIntegerAttribute("NM"))) {
                return o1.getIntegerAttribute("NM") < o2.getIntegerAttribute("NM") ? -1 : 1;
            }

            if (!o1.getCigar().equals(o2.getCigar())) {
                return o1.getCigar().numCigarElements() < o2.getCigar().numCigarElements() ? -1 : 1;
            }

            return 0;
        });

        return srs;
    }
}
