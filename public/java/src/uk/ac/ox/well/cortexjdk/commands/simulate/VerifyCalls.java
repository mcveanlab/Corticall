package uk.ac.ox.well.cortexjdk.commands.simulate;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.sw.SmithWaterman;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableWriter;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class VerifyCalls extends Module {
    @Argument(fullName="kmerTable", shortName="kt", doc="Kmer table")
    public File KMER_TABLE;

    @Argument(fullName = "variantTable", shortName = "vt", doc = "Variant table", required = false)
    public File VARIANT_TABLE;

    @Argument(fullName="rois", shortName="r", doc="ROIs")
    public CortexGraph ROI;

    @Argument(fullName="vcf", shortName="v", doc="VCF of novel variants")
    public ArrayList<VCFFileReader> VCFS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CanonicalKmer, Map<String, String>> kmerMap = new HashMap<>();
        Map<CanonicalKmer, Set<Integer>> kmerIds = new HashMap<>();
        Map<Integer, Map<String, String>> variantMap = new TreeMap<>();

        TableReader ktr = new TableReader(KMER_TABLE, "id", "length", "kmerIndex", "kmer");
        for (Map<String, String> te : ktr) {
            CanonicalKmer ck = new CanonicalKmer(te.get("kmer"));
            kmerMap.put(ck, te);

            if (!kmerIds.containsKey(ck)) {
                kmerIds.put(ck, new TreeSet<>());
            }

            kmerIds.get(ck).add(Integer.valueOf(te.get("id")));
        }

        if (VARIANT_TABLE != null) {
            TableReader vtr = new TableReader(VARIANT_TABLE);
            for (Map<String, String> te : vtr) {
                int variantId = Integer.valueOf(te.get("index"));
                if (variantId > -1) {
                    variantMap.put(variantId, te);
                }
            }
        }

        Map<Integer, VariantContext> matchedVCs = new HashMap<>();

        for (VCFFileReader vcf : VCFS) {
            for (Integer vid : variantMap.keySet()) {
                if (!variantMap.get(vid).get("type").contains("NAHR")) {
                    log.info("{} {}", vid, variantMap.get(vid));

                    double acc = 0.0;
                    VariantContext best = null;
                    String knownHap = (variantMap.get(vid).get("sleft") + variantMap.get(vid).get("new") + variantMap.get(vid).get("sright")).toUpperCase();

                    for (VariantContext vc : vcf) {
                        SmithWaterman sw = new SmithWaterman();
                        String[] a = sw.getAlignment(vc.getAttributeAsString("CHILD_HAP", ""), knownHap);

                        double numBases = 0.0;
                        double numMatches = 0.0;

                        int leftStart = 0;
                        for (int i = 0; i < a[0].length() && a[0].charAt(i) == '-'; i++) {
                            leftStart++;
                        }

                        int rightStart = a[0].length();
                        for (int i = a[0].length() - 1; i >= 0 && a[0].charAt(i) == '-'; i--) {
                            rightStart--;
                        }

                        for (int i = leftStart; i < rightStart; i++) {
                            numBases++;
                            if (a[0].charAt(i) == a[1].charAt(i)) {
                                numMatches++;
                            }
                        }

                        /*
                        log.info("{} {}", leftStart, rightStart);
                        log.info("{}", a[0]);
                        log.info("{}", a[1]);
                        log.info("{} {}", numMatches / numBases, acc);
                        */

                        if (numMatches / numBases > acc) {
                            acc = numMatches / numBases;
                            best = vc;
                        }
                    }

                    if (best != null) {
                        SmithWaterman sw = new SmithWaterman();
                        String[] a = sw.getAlignment(best.getAttributeAsString("CHILD_HAP", ""), knownHap);

                        log.info("{} {}", vid, acc);
                        log.info("{}", a[0]);
                        log.info("{}", a[1]);
                        log.info("{}", best);

                        matchedVCs.put(vid, best);
                    }
                }
            }

            for (VariantContext vc : vcf) {
                if (!matchedVCs.values().contains(vc)) {
                    log.info("{}", vc);
                }
            }
        }

        /*
        Set<CanonicalKmer> used = new HashSet<>();
        for (CortexRecord cr : ROI) {
            used.add(cr.getCanonicalKmer());
        }

        int kmerSize = kmerMap.keySet().iterator().next().length();

        Map<Integer, Map<String, String>> stats = new TreeMap<>();

        ReferenceSequence rseq;
        while ((rseq = FASTA.nextSequence()) != null) {
            String seq = rseq.getBaseString();

            Set<Integer> ids = new HashSet<>();
            Set<CanonicalKmer> kmers = new HashSet<>();
            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                CanonicalKmer ck = new CanonicalKmer(seq.substring(i, i + kmerSize));
                kmers.add(ck);

                if (kmerMap.containsKey(ck)) {
                    ids.add(Integer.valueOf(kmerMap.get(ck).get("id")));
                }
            }

            if (ids.size() > 0) {
                for (int variantId : ids) {
                    int numExp = 0, numFound = 0;
                    if (variantMap.containsKey(variantId)) {
                        Map<String, String> vm = variantMap.get(variantId);

                        String contigWithNewAllele = (vm.get("sleft") + (vm.get("new").equals(".") ? "" : vm.get("new")) + vm.get("sright")).toUpperCase();
                        for (int i = 0; i <= contigWithNewAllele.length() - ROI.getKmerSize(); i++) {
                            String sk = contigWithNewAllele.substring(i, i + ROI.getKmerSize());
                            CanonicalKmer ak = new CanonicalKmer(sk);

                            if (used.contains(ak)) {
                                numExp++;
                                if (kmers.contains(ak)) {
                                    numFound++;
                                }
                            }
                        }

                        String type = vm.get("type");
                        int alleleLength = vm.get("new").equals(".") ? 0 : vm.get("new").length();

                        //out.println(Joiner.on("\t").join(rseq.getName().split(" ")[0], type, alleleLength, numFound, numExp, numFound != numExp ? "incomplete" : "complete"));

                        Map<String, String> entry = new LinkedHashMap<>();
                        entry.put("variantId", String.valueOf(variantId));
                        entry.put("type", type);
                        entry.put("alleleLength", String.valueOf(alleleLength));
                        entry.put("numFound", String.valueOf(numFound));
                        entry.put("numExp", String.valueOf(numExp));
                        entry.put("completeness", numFound != numExp ? "incomplete" : "complete");

                        if (!stats.containsKey(variantId) || Integer.valueOf(stats.get(variantId).get("numFound")) < numFound) {
                            stats.put(variantId, entry);
                        }
                    }
                }
            }
        }

        TableWriter tw = new TableWriter(out);
        for (Map<String, String> te : stats.values()) {
            tw.addEntry(te);
        }
        */
    }
}
