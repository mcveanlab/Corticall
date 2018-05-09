package uk.ac.ox.well.cortexjdk.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableWriter;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class VerifyAlleleCompleteness extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public FastaSequenceFile FASTA;

    @Argument(fullName="kmerTable", shortName="kt", doc="Kmer table")
    public File KMER_TABLE;

    @Argument(fullName = "variantTable", shortName = "vt", doc = "Variant table", required = false)
    public File VARIANT_TABLE;

    @Argument(fullName="rois", shortName="r", doc="ROIs")
    public CortexGraph ROI;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CanonicalKmer, Map<String, String>> kmerMap = new HashMap<>();
        Map<CanonicalKmer, Set<Integer>> kmerIds = new HashMap<>();
        Map<Integer, Map<String, String>> variantMap = new HashMap<>();

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
    }
}
