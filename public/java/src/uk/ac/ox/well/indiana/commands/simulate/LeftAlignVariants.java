package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.util.HashSet;
import java.util.Set;

public class LeftAlignVariants extends Module {
    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Argument(fullName="reference", shortName="r", doc="Reference")
    public IndexedFastaSequenceFile REF;

    @Output
    public File out;

    private VariantContext leftShift(VariantContext vc) {
        if (vc.isIndel() || vc.isMNP()) {
            String r = vc.getReference().getBaseString();
            String a = vc.getAlternateAllele(0).getBaseString();
            String prevBase = new String(REF.getSubsequenceAt(vc.getChr(), vc.getStart(), vc.getStart()).getBases());
            String nextBase = (a.length() > r.length()) ? String.valueOf(a.charAt(a.length() - 1)) : String.valueOf(r.charAt(r.length() - 1));

            while (prevBase.equals(nextBase)) {
                String newr = new String(REF.getSubsequenceAt(vc.getChr(), vc.getStart(), vc.getStart() + r.length()).getBases());
                String newa = prevBase + ((a.length() - 1 < 0) ? a.substring(0, a.length() - 1) : "");

                if (newr.equals(newa)) {
                    break;
                } else {
                    vc = (new VariantContextBuilder(vc))
                            .start(vc.getStart() - 1)
                            .stop(vc.getStart() + (vc.getEnd() - vc.getStart()))
                            .alleles(newr, newa)
                            .attribute("shifted", true)
                            .make();

                    r = vc.getReference().getBaseString();
                    a = vc.getAlternateAllele(0).getBaseString();
                    prevBase = new String(REF.getSubsequenceAt(vc.getChr(), vc.getStart() - 1, vc.getStart()).getBases());
                    nextBase = (a.length() > r.length()) ? String.valueOf(a.charAt(a.length() - 1)) : String.valueOf(r.charAt(r.length() - 1));
                }
            }
        }

        return vc;
    }
    @Override
    public void execute() {
        VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder();
        vcwb.setOutputFile(out);
        vcwb.unsetOption(Options.INDEX_ON_THE_FLY);
        vcwb.setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        VariantContextWriter vcw = vcwb.build();

        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        VCFHeader header = new VCFHeader(headerLines);
        header.setSequenceDictionary(REF.getSequenceDictionary());

        vcw.writeHeader(header);

        log.info("Shifting variants...");
        int numShifted = 0, numTotal = 0;
        for (VariantContext vc : VCF) {
            if (vc.getStart() >= 0 && vc.getEnd() < REF.getSequence(vc.getChr()).length()) {
                vc = leftShift(vc);
            }
            vcw.add(vc);

            if (vc.hasAttribute("shifted")) {
                numShifted++;
            }
            numTotal++;
        }
        log.info("  {}/{} variants shifted", numShifted, numTotal);

        vcw.close();
    }
}
