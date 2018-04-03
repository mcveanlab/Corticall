package uk.ac.ox.well.cortexjdk.commands.discover.eval;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.containers.ContainerUtils;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableWriter;
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

    @Output
    public File out;

    @Override
    public void execute() {
        Map<CanonicalKmer, Set<VariantContext>> kmerMap = new HashMap<>();
        Map<VariantContext, Set<VariantContext>> recovered = new LinkedHashMap<>();

        for (VariantContext kvc : VCF_KNOWN) {
            recovered.put(kvc, new HashSet<>());

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

        Set<VariantContext> unused = new HashSet<>();
        for (VariantContext nvc : VCF_NOVEL) {
            Map<VariantContext, Integer> knownVcs = new HashMap<>();

            String childHap = nvc.getAttributeAsString("CHILD_HAP", "");
            for (int i = 0; i <= childHap.length() - KMER_SIZE; i++) {
                String sk = childHap.substring(i, i + KMER_SIZE);
                CanonicalKmer ck = new CanonicalKmer(sk);

                if (kmerMap.containsKey(ck)) {
                    //knownVcs.addAll(kmerMap.get(ck));

                    for (VariantContext vc : kmerMap.get(ck)) {
                        ContainerUtils.increment(knownVcs, vc);
                    }
                }
            }

            VariantContext bestVc = null;
            int bestCount = 0;
            for (VariantContext knownVc : knownVcs.keySet()) {
                if (knownVcs.get(knownVc) > bestCount) {
                    bestVc = knownVc;
                    bestCount = knownVcs.get(knownVc);
                }
            }

            if (bestVc != null) {
                recovered.get(bestVc).add(nvc);
            } else {
                unused.add(nvc);
            }

            log.info("call : {}", nvc);
            log.info("known: {} {}", bestCount, bestVc);
            log.info("");
        }

        TableWriter tw = new TableWriter(out);

        Map<String, Pair<Integer, Integer>> typeCount = new HashMap<>();
        for (VariantContext vc : recovered.keySet()) {
            Map<String, String> te = new LinkedHashMap<>();
            te.put("KNOWN_INDEX", vc.getAttributeAsString("INDEX", "-1"));
            te.put("KNOWN_TYPE", vc.getAttributeAsString("SIM_TYPE", "UNKNOWN"));
            te.put("KNOWN_CHR", vc.getContig());
            te.put("KNOWN_START", String.valueOf(vc.getStart()));
            te.put("KNOWN_END", String.valueOf(vc.getEnd()));
            te.put("KNOWN_LENGTH", String.valueOf(Math.max(vc.getReference().length(), vc.getAltAlleleWithHighestAlleleCount().length())));

            String type = vc.getAttributeAsString("SIM_TYPE", "UNKNOWN");

            if (!typeCount.containsKey(type)) {
                typeCount.put(type, Pair.create(0, 0));
            }

            int knownCount = typeCount.get(type).getFirst() + 1;
            int novelCount = typeCount.get(type).getSecond();

            if (recovered.get(vc).size() > 0) {
                Interval kit = new Interval(vc.getContig(), vc.getStart() - 10, vc.getEnd() + 10);

                for (VariantContext nvc : recovered.get(vc)) {
                    Interval nit = new Interval(nvc.getContig(), nvc.getStart() - 10, nvc.getEnd() + 10);

                    Map<String, String> te2 = new LinkedHashMap<>(te);
                    te2.put("NOVEL_PIECES", String.valueOf(recovered.get(vc).size()));
                    te2.put("NOVEL_TYPE", nvc.getType().name());
                    te2.put("NOVEL_CHR", nvc.getContig());
                    te2.put("NOVEL_START", String.valueOf(nvc.getStart()));
                    te2.put("NOVEL_END", String.valueOf(nvc.getEnd()));
                    te2.put("NOVEL_COMPLETE", (nvc.getAltAlleleWithHighestAlleleCount().getDisplayString().startsWith(".") || nvc.getAltAlleleWithHighestAlleleCount().getDisplayString().endsWith(".")) ? "FALSE" : "TRUE");
                    te2.put("NOVEL_COORDS_CORRECT", kit.getIntersectionLength(nit) > 0 ? "TRUE" : "FALSE");

                    tw.addEntry(te2);
                }

                novelCount++;
            } else {
                Map<String, String> te2 = new LinkedHashMap<>(te);
                te2.put("NOVEL_PIECES", "NA");
                te2.put("NOVEL_TYPE", "NA");
                te2.put("NOVEL_CHR", "NA");
                te2.put("NOVEL_START", "NA");
                te2.put("NOVEL_END", "NA");
                te2.put("NOVEL_COMPLETE", "NA");
                te2.put("NOVEL_COORDS_CORRECT", "NA");

                tw.addEntry(te2);
            }

            typeCount.put(type, Pair.create(knownCount, novelCount));
        }

        for (VariantContext nvc : unused) {
            Map<String, String> te = new LinkedHashMap<>();
            te.put("KNOWN_INDEX", "NA");
            te.put("KNOWN_TYPE", "NA");
            te.put("KNOWN_CHR", "NA");
            te.put("KNOWN_START", "NA");
            te.put("KNOWN_END", "NA");
            te.put("KNOWN_LENGTH", "NA");
            te.put("NOVEL_PIECES", "1");
            te.put("NOVEL_TYPE", nvc.getType().name());
            te.put("NOVEL_CHR", nvc.getContig());
            te.put("NOVEL_START", String.valueOf(nvc.getStart()));
            te.put("NOVEL_END", String.valueOf(nvc.getEnd()));
            te.put("NOVEL_COMPLETE", (nvc.getAltAlleleWithHighestAlleleCount().getDisplayString().startsWith(".") || nvc.getAltAlleleWithHighestAlleleCount().getDisplayString().endsWith(".")) ? "FALSE" : "TRUE");
            te.put("NOVEL_COORDS_CORRECT", "NA");

            tw.addEntry(te);
        }

        for (String type : typeCount.keySet()) {
            log.info("{} {} {}", type, typeCount.get(type).getFirst(), typeCount.get(type).getSecond());
        }
    }
}
