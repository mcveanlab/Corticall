package uk.ac.ox.well.cortexjdk.playground.eval;

import com.google.common.base.Joiner;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.commands.simulate.generators.GeneratedVariant;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class CompareCalls extends Module {
    @Argument(fullName="truth", shortName="t", doc="Truth table")
    public ArrayList<File> TRUTH_TABLES;

    @Argument(fullName="calls", shortName="c", doc="Calls")
    public ArrayList<File> CALLS;

    @Argument(fullName="references", shortName="R", doc="Reference(s)", required=false)
    public HashMap<String, IndexedReference> REFERENCES;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, Map<String, Integer>> stats = createStatsTable();

        Set<GeneratedVariant> allGvs = new HashSet<>();
        Set<VariantContext> allCvs = new HashSet<>();

        for (int i = 0; i < TRUTH_TABLES.size(); i++) {
            Map<GeneratedVariant, Map<VariantContext, Double>> m = compareHaplotypes(TRUTH_TABLES.get(i), CALLS.get(i));

            Set<VariantContext> all = new HashSet<>();
            for (GeneratedVariant gv : m.keySet()) {
                all.addAll(m.get(gv).keySet());
            }

            Set<VariantContext> usedV = new HashSet<>();
            Set<GeneratedVariant> usedG = new HashSet<>();

            allGvs.addAll(m.keySet());
            for (GeneratedVariant gv : m.keySet()) {
                double bestScore = 0.0;
                VariantContext bestCv = null;

                for (VariantContext cv : m.get(gv).keySet()) {
                    if (!usedV.contains(cv)) {
                        double score = m.get(gv).get(cv);
                        if (score >= bestScore) {
                            bestScore = score;
                            bestCv = cv;
                        }
                    }
                }

                double secondBestScore = 0.0;
                VariantContext secondBestCv = null;

                for (VariantContext cv : m.get(gv).keySet()) {
                    if (!usedV.contains(cv)) {
                        double score = m.get(gv).get(cv);
                        if (score >= secondBestScore && score < bestScore) {
                            secondBestScore = score;
                            secondBestCv = cv;
                        }
                    }
                }

                double totalScore = bestScore;

                if (bestCv != null && bestCv.isSymbolic() && bestScore < 0.5 && secondBestScore > 0.4) {
                    totalScore += secondBestScore;
                    usedV.add(secondBestCv);
                }

                if (totalScore > 0.8) {
                    int length = gv.newAllele.length() == gv.oldAllele.length() ? gv.newAllele.length() : Math.abs(gv.newAllele.length() - gv.oldAllele.length());

                    if (gv.getType().equals("SNV")) {
                        increment(stats, "SNV", 0, "N");
                        if (bestCv != null) {
                            increment(stats, "SNV", 0, "TP");
                        } else {
                            increment(stats, "SNV", 0, "FN");
                        }
                    } else if (gv.getType().equals("STR_EXP")) {
                        increment(stats, "STR_EXP", 0, "N");
                        if (bestCv != null) {
                            increment(stats, "STR_EXP", 0, "TP");
                        } else {
                            increment(stats, "STR_EXP", 0, "FN");
                        }
                    } else if (gv.getType().equals("STR_CON")) {
                        increment(stats, "STR_CON", 0, "N");
                        if (bestCv != null) {
                            increment(stats, "STR_CON", 0, "TP");
                        } else {
                            increment(stats, "STR_CON", 0, "FN");
                        }
                    } else if (gv.getType().equals("NAHR-INS")) {
                        increment(stats, "NAHR-INS", 0, "N");
                        if (bestCv != null) {
                            increment(stats, "NAHR-INS", 0, "TP");
                        } else {
                            increment(stats, "NAHR-INS", 0, "FN");
                        }
                    } else if (gv.getType().equals("TD")) {
                        increment(stats, "TD", length, "N");
                        if (bestCv != null) {
                            increment(stats, "TD", length, "TP");
                        } else {
                            increment(stats, "TD", length, "FN");
                        }
                    } else if (gv.getType().equals("INV")) {
                        increment(stats, "INV", length, "N");
                        if (bestCv != null) {
                            increment(stats, "INV", length, "TP");
                        } else {
                            increment(stats, "INV", length, "FN");
                        }
                    } else if (gv.getType().equals("DEL") || gv.getType().equals("NAHR-DEL")) {
                        increment(stats, "DEL", length, "N");
                        if (bestCv != null) {
                            increment(stats, "DEL", length, "TP");
                        } else {
                            increment(stats, "DEL", length, "FN");
                        }
                    } else if (gv.getType().equals("MNP")) {
                        increment(stats, "MNP", length, "N");
                        if (bestCv != null) {
                            increment(stats, "MNP", length, "TP");
                        } else {
                            increment(stats, "MNP", length, "FN");
                        }
                    } else if (gv.getType().equals("INS")) {
                        increment(stats, "INS", length, "N");
                        if (bestCv != null) {
                            increment(stats, "INS", length, "TP");
                        } else {
                            increment(stats, "INS", length, "FN");
                        }
                    } else {
                        log.info("");
                    }

                    if (bestCv != null) {
                        usedV.add(bestCv);
                        usedG.add(gv);
                    }
                }
            }

            for (VariantContext cv : all) {
                if (!usedV.contains(cv)) {
                    double bestScore = 0.0;
                    GeneratedVariant bestGv = null;

                    for (GeneratedVariant gv : m.keySet()) {
                        if (usedG.contains(gv)) {
                            double score = computeSimilarity(gv, cv);
                            if (score > bestScore) {
                                bestScore = score;
                                bestGv = gv;
                            }
                        }
                    }

//                    if (bestGv == null || bestScore < 0.5) {
//                        log.info("{} {} {}", bestScore, bestGv, cv);
//                        int k = 0;
//                    }
                }
            }

//            for (GeneratedVariant gv : m.keySet()) {
//                for (VariantContext cv : m.get(gv).keySet()) {
//                    if (!used.contains(cv)) {
//                        int length = cv.isMNP() ? cv.getReference().getBaseString().length() : Math.abs(cv.getAlternateAllele(0).length() - cv.getReference().getBaseString().length());
//
//                        if (cv.isSNP()) { increment(stats, "SNV", 0, "FP"); }
//                        else if (cv.isSimpleInsertion()) { increment(stats, "INS", length, "FP"); }
//                        else if (cv.isMNP()) {
//                            if (cv.getReference().getBaseString().equals(SequenceUtils.reverseComplement(cv.getAlternateAllele(0).getBaseString()))) {
//                                increment(stats, "INV", length, "FP");
//                            } else {
//                                increment(stats, "MNP", length, "FP");
//                            }
//                        }
//                        else if (cv.getAlternateAllele(0).getDisplayString().contains(":") && !cv.getAlternateAllele(0).getDisplayString().contains("TANDEM")) {
//                            if (cv.getAlternateAllele(0).getDisplayString().contains(cv.getContig())) {
//                                increment(stats, "UNKNOWN", 0, "FP");
//                            } else {
//                                increment(stats, "NAHR-INS", 0, "FP");
//                                Pair<String, Long> l = convertAltToLocus(cv);
//                                int j = 0;
//                            }
//                        } else {
//                            increment(stats, "UNKNOWN", 0, "FP");
//                        }
//
//                        used.add(cv);
//                    }
//                }
//            }
        }

        log.info("N: {} {}", allGvs.size(), allCvs.size());

        printTable(stats);
    }

    private void printTable(Map<String, Map<String, Integer>> stats) {
        out.println(Joiner.on("\t").join(
                String.format("%-20s", "type"),
                String.format("%-7s", "N"),
                String.format("%-7s", "TP"),
                String.format("%-7s", "FN"),
                String.format("%-7s", "FP"),
                String.format("%-7s", "F1")
        ));

        int expVariants = 0, allVariants = 0;

        for (String key : stats.keySet()) {
            double tp = (double) stats.get(key).get("TP");
            double fn = (double) stats.get(key).get("FN");
            double fp = (double) stats.get(key).get("FP");

            expVariants += tp + fn;
            allVariants += tp + fp;

            double f1 = (tp == 0 && fp == 0 && fn == 0) ? 0.0 : (2*tp) / ((2*tp) + fp + fn);
            out.println(Joiner.on("\t").join(
                    String.format("%-20s", key.replaceAll(" ", "_")),
                    String.format("%-7d", stats.get(key).get("N")),
                    String.format("%-7d", stats.get(key).get("TP")),
                    String.format("%-7d", stats.get(key).get("FN")),
                    String.format("%-7d", stats.get(key).get("FP")),
                    String.format("%-7s", String.format("%.2f", f1))
            ));
        }

        log.info("Total: {} {}", expVariants - 1, allVariants);
    }

    private void increment(Map<String, Map<String, Integer>> stats, String type, int length, String key) {
        String L = "";
        if (length == 0) { L = ""; }
        else if (length <= 100) { L = " (1-100)"; }
        else if (length <= 500) { L = " (101-500)"; }
        else { L = " (501-1000)"; }

        stats.get(type + L).put(key, stats.get(type + L).get(key) + 1);
    }

    private Map<String, Map<String, Integer>> createStatsTable() {
        Map<String, Map<String, Integer>> stats = new LinkedHashMap<>();

        for (String key : Arrays.asList("UNKNOWN", "SNV", "STR_EXP", "STR_CON", "NAHR-INS")) {
            Map<String, Integer> h = new LinkedHashMap<>();
            h.put("N", 0);
            h.put("TP", 0);
            h.put("FP", 0);
            h.put("FN", 0);

            stats.put(key, h);
        }

        for (String key : Arrays.asList("TD", "INV", "DEL", "MNP", "INS")) {
            for (String l : Arrays.asList("(1-100)", "(101-500)", "(501-1000)")) {
                Map<String, Integer> h = new LinkedHashMap<>();
                h.put("N", 0);
                h.put("TP", 0);
                h.put("FP", 0);
                h.put("FN", 0);

                stats.put(key + " " + l, h);
            }
        }

        return stats;
    }

    private Map<GeneratedVariant, Map<VariantContext, Double>> compareHaplotypes(File truthTableFile, File variantFile) {
        Set<VariantContext> cvs = loadCalledVariants(variantFile);
        Set<GeneratedVariant> gvs = loadSimulatedVariants(truthTableFile);

        Map<GeneratedVariant, Map<VariantContext, Double>> m = new HashMap<>();

        for (GeneratedVariant gv : gvs) {
            m.put(gv, new HashMap<>());
            for (VariantContext cv : cvs) {
                m.get(gv).put(cv, computeSimilarity(gv, cv));
            }
            int i = 0;
        }

        return m;
    }

    private Pair<String, Long> convertAltToLocus(VariantContext cv) {
        String q = cv.getAlternateAllele(0).getDisplayString()
                .replaceAll("(HB3|DD2):", "")
                .replaceAll("[\\[\\]]", "")
                .replaceAll("^[ACGTN]+", "")
                .replaceAll("[ACGTN]+$", "")
                .replaceAll("-.*", "");

        String[] qa = q.split(":");
        String chr2 = qa[0];
        long start2 = Integer.parseInt(qa[1]);

        return new Pair<>(chr2, start2);
    }

    private double computeSimilarity(GeneratedVariant gv, VariantContext cv) {
        double score = 0.0;
        if (cv.getAlternateAllele(0).getDisplayString().contains(":") && !cv.getAlternateAllele(0).getDisplayString().contains("TANDEM")) {
            Pair<String, Long> l = convertAltToLocus(cv);
            String chr2 = l.getFirst();
            long start2 = l.getSecond();

            int nahrSeqIndexGV1 = gv.start.getContig().contains("HB3") ?
                    REFERENCES.get("HB3").getReferenceSequence().getSequence(gv.start.getContig()).getContigIndex() :
                    REFERENCES.get("DD2").getReferenceSequence().getSequence(gv.start.getContig()).getContigIndex();

            int nahrSeqIndexGV2 = gv.stop.getContig().contains("HB3") ?
                    REFERENCES.get("HB3").getReferenceSequence().getSequence(gv.stop.getContig()).getContigIndex() :
                    REFERENCES.get("DD2").getReferenceSequence().getSequence(gv.stop.getContig()).getContigIndex();

            int nahrSeqIndexV1 = cv.getContig().contains("HB3") ?
                    REFERENCES.get("HB3").getReferenceSequence().getSequence(cv.getContig()).getContigIndex() :
                    REFERENCES.get("DD2").getReferenceSequence().getSequence(cv.getContig()).getContigIndex();

            int nahrSeqIndexV2 = chr2.contains("HB3") ?
                    REFERENCES.get("HB3").getReferenceSequence().getSequence(chr2).getContigIndex() :
                    REFERENCES.get("DD2").getReferenceSequence().getSequence(chr2).getContigIndex();

            boolean sameSeqIndex = (nahrSeqIndexGV1 == nahrSeqIndexV1 && nahrSeqIndexGV2 == nahrSeqIndexV2) || (nahrSeqIndexGV1 == nahrSeqIndexV2 && nahrSeqIndexGV2 == nahrSeqIndexV1);

            double closestLeft1 = Math.abs(gv.start.getStart() - cv.getStart());
            double closestLeft2 = Math.abs(gv.start.getStart() - start2);

            double closestRight1 = Math.abs(gv.stop.getStart() - cv.getStart());
            double closestRight2 = Math.abs(gv.stop.getStart() - start2);

            double closestLeft = Math.min(closestLeft1, closestLeft2);
            double closestRight = Math.min(closestRight1, closestRight2);

            score += (nahrSeqIndexV1 == nahrSeqIndexGV1 || nahrSeqIndexV1 == nahrSeqIndexGV2) ? (5000.0 - closestLeft)/10000.0 : 0.0;
            score += (nahrSeqIndexV2 == nahrSeqIndexGV1 || nahrSeqIndexV2 == nahrSeqIndexGV1) ? (5000.0 - closestRight)/10000.0 : 0.0;

            //return (10000.0 - (closestLeft + closestRight)) / 10000.0;
        } else {
            String simulatedHaplotype = (gv.seedLeft + gv.newAllele + gv.seedRight).replaceAll("\\.", "").toUpperCase();
            String calledHaplotype = "";

            IndexedReference ref = cv.getContig().contains("HB3") ? REFERENCES.get("HB3") : REFERENCES.get("DD2");

            String hleft = ref.getReferenceSequence().getSubsequenceAt(cv.getContig(), cv.getStart() - 100, cv.getStart() - 1).getBaseString();
            String hright = ref.getReferenceSequence().getSubsequenceAt(cv.getContig(), cv.getEnd() + 1, cv.getEnd() + 100).getBaseString();

            calledHaplotype = hleft + cv.getAlternateAllele(0).getDisplayString() + hright;

            score = compareKmerSets(simulatedHaplotype, calledHaplotype, 11);
        }

        return score;
    }

    private double compareKmerSets(String simulatedHaplotype, String calledHaplotype, int kmerSize) {
        Set<CanonicalKmer> kmersG = kmerize(simulatedHaplotype, kmerSize);
        Set<CanonicalKmer> kmersC = kmerize(calledHaplotype, kmerSize);
        Set<CanonicalKmer> kmersA = new HashSet<>();
        kmersA.addAll(kmersG);
        kmersA.addAll(kmersC);

        int numKmersOverlap = 0, numKmers = 0;
        for (CanonicalKmer kmer : kmersA) {
            if (kmersG.contains(kmer)) {
                numKmers++;

                if (kmersC.contains(kmer)) {
                    numKmersOverlap++;
                }
            }
        }

        double pctOverlap = (double) numKmersOverlap / (double) numKmers;

        return pctOverlap;
    }

    private Set<CanonicalKmer> kmerize(String seq, int kmerSize) {
        Set<CanonicalKmer> kmers = new HashSet<>();

        for (int i = 0; i <= seq.length() - kmerSize; i++) {
            String sk = seq.substring(i, i + kmerSize).toUpperCase();

            kmers.add(new CanonicalKmer(sk));
        }

        return kmers;
    }

    private Set<GeneratedVariant> loadSimulatedVariants(File truthTableFile) {
        Set<GeneratedVariant> gvs = new HashSet<>();

        TableReader truthTr = new TableReader(truthTableFile);
        for (Map<String, String> tr : truthTr) {
            if (!tr.get("type").equals("RECOMB")) {
                Interval start, stop;
                if (tr.get("refStart").contains(":")) {
                    String[] refStartPieces = tr.get("refStart").split("[:-]");
                    start = new Interval(refStartPieces[0], Integer.valueOf(refStartPieces[1]), Integer.valueOf(refStartPieces[2]), refStartPieces.length == 3, ".");

                    String[] refStopPieces = tr.get("refStop").split("[:-]");
                    stop = new Interval(refStopPieces[0], Integer.valueOf(refStopPieces[1]), Integer.valueOf(refStopPieces[2]), refStopPieces.length == 3, ".");
                } else {
                    start = new Interval(tr.get("refChr"), Integer.valueOf(tr.get("refStart")), Integer.valueOf(tr.get("refStart")));
                    stop = new Interval(tr.get("refChr"), Integer.valueOf(tr.get("refStop")), Integer.valueOf(tr.get("refStop")));
                }

                GeneratedVariant gv = new GeneratedVariant(
                        tr.get("type"),
                        Integer.valueOf(tr.get("chr")) - 1,
                        Integer.valueOf(tr.get("start")),
                        tr.get("old"),
                        tr.get("new"),
                        tr.get("sleft"),
                        tr.get("sright"),
                        tr.get("parent"),
                        start,
                        stop,
                        Integer.valueOf(tr.get("index"))
                );

                gvs.add(gv);
            }
        }

        return gvs;
    }

    private Set<VariantContext> loadCalledVariants(File variantFile) {
        VCFFileReader vcf = new VCFFileReader(variantFile, false);
        Set<VariantContext> vcs = new HashSet<>();
        Set<String> ignore = new HashSet<>();
        for (VariantContext vc : vcf) {
            IndexedReference ref = (vc.getContig().contains("HB3")) ? REFERENCES.get("HB3") : REFERENCES.get("DD2");

            if (!vc.getContig().contains("contig") &&
                !vc.getAlternateAllele(0).getDisplayString().contains("contig") &&
                vc.getStart() + 100 < ref.getReferenceSequence().getSequence(vc.getContig()).length() &&
                !vc.isFiltered() &&
                !ignore.contains(vc.getID()) &&
                !ignore.contains(vc.getAttributeAsString("MATEID", String.valueOf(vc.hashCode())))
            ) {
                vcs.add(vc);

                if (vc.hasAttribute("MATEID")) {
                    ignore.add(vc.getAttributeAsString("MATEID", String.valueOf(vc.hashCode())));
                    ignore.add(vc.getID());
                }
            }
        }

        return vcs;
    }
}

