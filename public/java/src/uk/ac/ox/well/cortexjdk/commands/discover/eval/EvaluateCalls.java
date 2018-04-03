package uk.ac.ox.well.cortexjdk.commands.discover.eval;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;

import java.io.File;
import java.util.*;

public class EvaluateCalls extends Module {
    @Argument(fullName="known", shortName="k", doc="VCF of known variants")
    public VCFFileReader VCF_KNOWN;

    @Argument(fullName="novel", shortName="n", doc="VCF of novel variants")
    public VCFFileReader VCF_NOVEL;

    @Argument(fullName="contigs", shortName="c", doc="Contigs", required=false)
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="background", shortName="b", doc="Background", required=false)
    public HashMap<String, IndexedReference> BACKGROUNDS;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Override
    public void execute() {
        Map<CanonicalKmer, Set<VariantContext>> kmerMap = new HashMap<>();

        for (VariantContext kvc : VCF_KNOWN) {
            String newHap = kvc.getAttributeAsString("NEW_HAP", "");
            for (int i = 0; i <= newHap.length() - KMER_SIZE; i++) {
                String sk = newHap.substring(i, i + KMER_SIZE);
                CanonicalKmer ck = new CanonicalKmer(sk);

                if (!kmerMap.containsKey(ck)) {
                    kmerMap.put(ck, new HashSet<>());
                }

                kmerMap.get(ck).add(kvc);
            }
        }

        for (VariantContext nvc : VCF_NOVEL) {
            Set<VariantContext> knownVcs = new HashSet<>();

            String childHap = nvc.getAttributeAsString("CHILD_HAP", "");
            for (int i = 0; i <= childHap.length() - KMER_SIZE; i++) {
                String sk = childHap.substring(i, i + KMER_SIZE);
                CanonicalKmer ck = new CanonicalKmer(sk);

                if (kmerMap.containsKey(ck)) {
                    knownVcs.addAll(kmerMap.get(ck));
                }
            }

            log.info("call : {}", nvc);
            for (VariantContext knownVc : knownVcs) {
                log.info("known: {}", knownVc);
            }
            log.info("");
        }
    }
}
