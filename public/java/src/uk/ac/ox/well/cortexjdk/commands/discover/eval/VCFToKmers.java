package uk.ac.ox.well.cortexjdk.commands.discover.eval;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;

import java.io.PrintStream;

public class VCFToKmers extends Module {
    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Argument(fullName="reference", shortName="R", doc="Reference")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 63;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        for (VariantContext vc : VCF) {
            String seqTotal = REF.getSubsequenceAt(vc.getContig(), vc.getStart() - KMER_SIZE, vc.getEnd() + KMER_SIZE).getBaseString();
            String seqBefore = REF.getSubsequenceAt(vc.getContig(), vc.getStart() - KMER_SIZE, vc.getStart() - 1).getBaseString();
            String seqAfter = REF.getSubsequenceAt(vc.getContig(), vc.getEnd() + 1, vc.getEnd() + KMER_SIZE).getBaseString();

            String seqOld = seqBefore + vc.getReference().getBaseString() + seqAfter;
            String seqNew = seqBefore + vc.getAltAlleleWithHighestAlleleCount().getBaseString() + seqAfter;

            /*
            log.info("{} {}", vc.getReference().getBaseString(), vc.getAltAlleleWithHighestAlleleCount().getBaseString());
            log.info("total: {}", seqTotal);
            log.info("  old: {}", seqOld);
            log.info("  new: {}", seqNew);
            */

            for (int i = 0; i <= seqNew.length() - KMER_SIZE; i++) {
                String sk = seqNew.substring(i, i + KMER_SIZE);
                CanonicalKmer ck = new CanonicalKmer(sk);

                out.println(Joiner.on("\t").join(vc.getContig(), vc.getStart(), vc.getReference().getDisplayString(), vc.getAltAlleleWithHighestAlleleCount().getDisplayString(), i, sk, ck));
            }
        }
    }
}
