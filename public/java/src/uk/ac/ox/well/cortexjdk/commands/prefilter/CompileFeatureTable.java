package uk.ac.ox.well.cortexjdk.commands.prefilter;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableWriter;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

public class CompileFeatureTable extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="rois", shortName="r", doc="ROIs")
    public CortexGraph ROIS;

    @Argument(fullName="feature", shortName="f", doc="Feature file")
    public TreeMap<String, CortexGraph> FEATURES;

    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="roisTruth", shortName="rt", doc="ROI truth")
    public CortexGraph ROIS_TRUTH;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<CanonicalKmer> truth = new HashSet<>();
        for (CortexRecord cr : ROIS_TRUTH) {
            truth.add(cr.getCanonicalKmer());
        }

        Map<CanonicalKmer, Map<String, Object>> featureTable = new LinkedHashMap<>();

        for (CortexRecord cr : ROIS) {
            featureTable.put(cr.getCanonicalKmer(), new HashMap<>());
        }

        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            String name = rseq.getName().split(" ")[0];
            String seq = rseq.getBaseString();

            Map<String, Object> featureEntry = new LinkedHashMap<>();
            Set<CanonicalKmer> novelKmersInPartition = new HashSet<>();
            for (int i = 0; i <= seq.length() - ROIS.getKmerSize(); i++) {
                CanonicalKmer ck = new CanonicalKmer(seq.substring(i, i + ROIS.getKmerSize()));

                if (featureTable.containsKey(ck)) {
                    int distanceFromTerminus = Math.min(i, seq.length() - ROIS.getKmerSize() - i);

                    featureEntry.put("partitionName", name);
                    featureEntry.put("partitionLength", seq.length());
                    featureEntry.put("distanceFromTerminus", distanceFromTerminus);
                    featureEntry.put("compressionRatio", SequenceUtils.computeCompressionRatio(ck));
                    featureEntry.put("truth", truth.contains(ck) ? 1 : 0);

                    novelKmersInPartition.add(ck);
                }
            }

            for (CanonicalKmer ck : novelKmersInPartition) {
                if (novelKmersInPartition.size() > (int) featureTable.get(ck).getOrDefault("numNovelsInPartition", 0)) {
                    featureTable.get(ck).putAll(featureEntry);
                    featureTable.get(ck).put("numNovelsInPartition", novelKmersInPartition.size());
                }
            }
        }

        for (String feature : FEATURES.keySet()) {
            Set<CanonicalKmer> cks = new HashSet<>();
            for (CortexRecord cr : FEATURES.get(feature)) {
                cks.add(cr.getCanonicalKmer());
            }

            for (CanonicalKmer ck : featureTable.keySet()) {
                featureTable.get(ck).put(feature, cks.contains(ck) ? 1 : 0);
            }
        }

        TableWriter tw = new TableWriter(out);

        for (CanonicalKmer ck : featureTable.keySet()) {
            Map<String, String> te = new LinkedHashMap<>();

            te.put("ck",                   ck.getKmerAsString());
            te.put("partitionName",        String.valueOf(featureTable.get(ck).getOrDefault("partitionName", "unknown")));
            te.put("partitionLength",      String.valueOf(featureTable.get(ck).getOrDefault("partitionLength", 0)));
            te.put("numNovelsInPartition", String.valueOf(featureTable.get(ck).getOrDefault("numNovelsInPartition", 0)));
            te.put("distanceFromTerminus", String.valueOf(featureTable.get(ck).getOrDefault("distanceFromTerminus", 0)));
            te.put("compressionRatio",     String.valueOf(featureTable.get(ck).getOrDefault("compressionRatio", 1.0f)));

            for (String feature : FEATURES.keySet()) {
                te.put(feature,            String.valueOf(featureTable.get(ck).getOrDefault(feature, 0)));
            }

            te.put("truth",                String.valueOf(featureTable.get(ck).getOrDefault("truth", 0)));

            tw.addEntry(te);
        }
    }
}
