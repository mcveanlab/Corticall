package uk.ac.ox.well.cortexjdk.playground.eval;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.commands.simulate.generators.GeneratedVariant;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
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

    @Argument(fullName="requireCorrectParent", shortName="P", doc="Require same parent")
    public Boolean REQUIRE_CORRECT_PARENT = false;

    @Argument(fullName="excludeDuplicatesInSyntenicLoci", shortName="X", doc="Exclude duplicate variants found in syntenic loci")
    public Boolean EXCLUDE_DUPLICATES_IN_SYNTENIC_LOCI = true;

    @Argument(fullName="proximalMatchRange", shortName="r", doc="Consider variants within specified range to be matches")
    public Integer PROXIMAL_MATCH_RANGE = 0;

    @Output
    public PrintStream out;

//    @Output(fullName="featureOut", shortName="fo", doc="Features out")
//    public PrintStream fout;

    @Override
    public void execute() {
        Map<String, Map<String, Integer>> stats = createStatsTable();

        Set<GeneratedVariant> allGvs = new HashSet<>();
        Set<VariantContext> allCvs = new HashSet<>();
        Set<VariantContext> unusedCVs = new HashSet<>();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing callsets")
                .maxRecord(TRUTH_TABLES.size())
                .message("callsets processed")
                .make(log);

        for (int i = 0; i < TRUTH_TABLES.size(); i++) {
            pm.update();

            Set<GeneratedVariant> gvs = loadSimulatedVariants(TRUTH_TABLES.get(i));
            Set<VariantContext> cvs = loadCalledVariants(CALLS.get(i));

            allGvs.addAll(gvs);
            allCvs.addAll(cvs);

            Map<GeneratedVariant, Map<VariantContext, Double>> ma = compareConcreteVariants(gvs, cvs);
            Map<GeneratedVariant, Map<Pair<VariantContext, VariantContext>, Double>> mb = compareSymbolicVariants(gvs, cvs);

            Set<VariantContext> usedV = new HashSet<>();

            Map<Integer, Map<GeneratedVariant, Set<String>>> nahrsFound = new HashMap<>();

            for (GeneratedVariant gv : gvs) {
                if (gv.type.equals("NAHR-INS")) {
                    if (!nahrsFound.containsKey(gv.index)) {
                        nahrsFound.put(gv.index, new HashMap<>());
                    }
                    nahrsFound.get(gv.index).put(gv, new HashSet<>());
                }

                double bestScore = 0.0;
                VariantContext bestCv = null;

                for (VariantContext cv : ma.get(gv).keySet()) {
                    double score = ma.get(gv).get(cv);
                    if (score >= bestScore) {
                        bestScore = score;
                        bestCv = cv;
                    }
                }

                double bestPairScore = 0.0;
                Pair<VariantContext, VariantContext> bestPairCv = null;
                if (bestCv == null) {
                    for (Pair<VariantContext, VariantContext> p : mb.get(gv).keySet()) {
                        double score = mb.get(gv).get(p);
                        if (score >= bestPairScore) {
                            bestPairScore = score;
                            bestPairCv = p;
                        }
                    }
                }

                if (PROXIMAL_MATCH_RANGE > 0 && bestScore < 0.8 && bestPairScore < 0.8) {
                    VariantContext closestCalledVariantStart = closestCalledVariant(gv, cvs, true);
                    VariantContext closestCalledVariantStop = closestCalledVariant(gv, cvs, false);

                    if (closestCalledVariantStart != null && closestCalledVariantStop != null) {
                        if (closestCalledVariantStart.equals(closestCalledVariantStop)) {
                            bestCv = closestCalledVariantStart;
                            bestScore = 0.9;
                        } else {
                            bestPairCv = new Pair<>(closestCalledVariantStart, closestCalledVariantStop);
                            bestPairScore = 0.9;
                        }
                    }
                }

                int length = gv.newAllele.length() == gv.oldAllele.length() ? gv.newAllele.length() : Math.abs(gv.newAllele.length() - gv.oldAllele.length());

                switch (gv.getType()) {
                    case "SNV":
                        increment(stats, "SNV", 0, "N");
                        if (bestScore >= 0.8 && bestCv != null) {
                            increment(stats, "SNV", 0, "TP");
                            usedV.add(bestCv);
                        } else {
                            increment(stats, "SNV", 0, "FN");
                        }
                        break;
                    case "STR_EXP":
                        increment(stats, "STR_EXP", 0, "N");
                        if (bestScore >= 0.8 && bestCv != null) {
                            increment(stats, "STR_EXP", 0, "TP");
                            usedV.add(bestCv);
                        } else {
                            increment(stats, "STR_EXP", 0, "FN");
                        }
                        break;
                    case "STR_CON":
                        increment(stats, "STR_CON", 0, "N");
                        if (bestScore >= 0.8 && bestCv != null) {
                            increment(stats, "STR_CON", 0, "TP");
                            usedV.add(bestCv);
                        } else {
                            increment(stats, "STR_CON", 0, "FN");
                        }
                        break;
                    case "NAHR-INS":
                        increment(stats, "NAHR-INS", 0, "N");
                        if (bestPairScore >= 0.8 && bestPairCv != null) {
                            increment(stats, "NAHR-INS", 0, "TP");
                            usedV.add(bestPairCv.getFirst());
                            usedV.add(bestPairCv.getSecond());

                            nahrsFound.get(gv.index).get(gv).add(bestPairCv.getFirst().getAttributeAsString("PARTITION_NAME", bestPairCv.getFirst().getAlternateAllele(0).getDisplayString()));
                            nahrsFound.get(gv.index).get(gv).add(bestPairCv.getSecond().getAttributeAsString("PARTITION_NAME", bestPairCv.getSecond().getAlternateAllele(0).getDisplayString()));
                        } else {
                            increment(stats, "NAHR-INS", 0, "FN");
                        }
                        break;
                    case "TD":
                        increment(stats, "TD", length, "N");
                        if (bestScore >= 0.8 && bestCv != null) {
                            increment(stats, "TD", length, "TP");
                            usedV.add(bestCv);
                        } else if (bestPairScore >= 0.8 && bestPairCv != null) {
                            increment(stats, "TD", length, "TP");
                            usedV.add(bestPairCv.getFirst());
                            usedV.add(bestPairCv.getSecond());
                        } else {
                            increment(stats, "TD", length, "FN");
                        }
                        break;
                    case "INV":
                        increment(stats, "INV", length, "N");
                        if (bestScore >= 0.8 && bestCv != null) {
                            increment(stats, "INV", length, "TP");
                            usedV.add(bestCv);
                        } else if (bestPairScore >= 0.8 && bestPairCv != null) {
                            increment(stats, "INV", length, "TP");
                            usedV.add(bestPairCv.getFirst());
                            usedV.add(bestPairCv.getSecond());
                        } else {
                            increment(stats, "INV", length, "FN");
                        }
                        break;
                    case "DEL":
                    case "NAHR-DEL":
                        increment(stats, "DEL", length, "N");
                        if (bestScore >= 0.8 && bestCv != null) {
                            increment(stats, "DEL", length, "TP");
                            usedV.add(bestCv);
                        } else if (bestPairScore >= 0.8 && bestPairCv != null) {
                            increment(stats, "DEL", length, "TP");
                            usedV.add(bestPairCv.getFirst());
                            usedV.add(bestPairCv.getSecond());
                        } else {
                            increment(stats, "DEL", length, "FN");
                        }
                        break;
                    case "MNP":
                        increment(stats, "MNP", length, "N");
                        if (bestScore >= 0.8 && bestCv != null) {
                            increment(stats, "MNP", length, "TP");
                            usedV.add(bestCv);
                        } else if (bestPairScore >= 0.8 && bestPairCv != null) {
                            increment(stats, "MNP", length, "TP");
                            usedV.add(bestPairCv.getFirst());
                            usedV.add(bestPairCv.getSecond());
                        } else {
                            increment(stats, "MNP", length, "FN");
                        }
                        break;
                    case "INS":
                        increment(stats, "INS", length, "N");
                        if (bestScore >= 0.8 && bestCv != null) {
                            increment(stats, "INS", length, "TP");
                            usedV.add(bestCv);
                        } else if (bestPairScore >= 0.8 && bestPairCv != null) {
                            increment(stats, "INS", length, "TP");
                            usedV.add(bestPairCv.getFirst());
                            usedV.add(bestPairCv.getSecond());
                        } else {
                            increment(stats, "INS", length, "FN");
                        }
                        break;
                    default:
                        break;
                }
            }

            for (Integer index : nahrsFound.keySet()) {
                boolean allFound = true;
                Set<String> partitionNames = new HashSet<>();
                for (GeneratedVariant gv : nahrsFound.get(index).keySet()) {
                    allFound &= nahrsFound.get(index).get(gv).size() > 0;
                    partitionNames.addAll(nahrsFound.get(index).get(gv));
                }

                increment(stats, "NAHR-COMPLETE", 0, "N");
                if (allFound && partitionNames.size() == 1) {
                    increment(stats, "NAHR-COMPLETE", 0, "TP");
                } else {
                    increment(stats, "NAHR-COMPLETE", 0, "FN");
                }
            }

            Set<VariantContext> moreUsedV = new HashSet<>();
            for (VariantContext cv : cvs) {
                if (!usedV.contains(cv)) {
                    for (VariantContext u : usedV) {
                        if (cv.getContig().equals(u.getContig()) && Math.abs(cv.getStart() - u.getStart()) <= 100) {
                            moreUsedV.add(cv);
                        }
                    }
                }
            }

            usedV.addAll(moreUsedV);

            Set<VariantContext> unusedV = new HashSet<>();
            for (VariantContext cv : cvs) {
                if (!usedV.contains(cv)) {
                    unusedV.add(cv);
                }
            }

            Set<String> seenIds = new HashSet<>();
            for (VariantContext cv : unusedV) {
                String id1 = cv.getID();
                String id2 = cv.getAttributeAsString("MATEID", cv.getID());

                Interval a = new Interval(cv.getContig(), cv.getStart() - 500, cv.getEnd() + cv.getAlternateAllele(0).length() + 500);

                boolean overlap = false;
                for (GeneratedVariant gv : gvs) {
                    Interval b = a;

                    if (cv.isSymbolic() && cv.getAlternateAllele(0).getDisplayString().contains(":") && !cv.getAlternateAllele(0).getDisplayString().contains("TANDEM")) {
                        b = convertAltToLocus(cv);
                    }

                    if (a.intersects(gv.start) || a.intersects(gv.stop) || b.intersects(gv.start) || b.intersects(gv.stop)) {
                        overlap = true;
                    }
                }

                if (!seenIds.contains(id1) && !seenIds.contains(id2) && !overlap) {
                    int length = cv.isMNP() ? cv.getReference().getBaseString().length() : Math.abs(cv.getAlternateAllele(0).length() - cv.getReference().getBaseString().length());

                    if (cv.isSNP()) {
                        increment(stats, "SNV", 0, "FP");
                        usedV.add(cv);

//                        fout.println(Joiner.on("\t").join(
//                                "SNV",
//                                "FP",
//                                cv.getAttributeAsDouble("SOR", 0.0),
//                                cv.getAttributeAsDouble("QD", 0.0),
//                                cv.getAttributeAsDouble("MQ", 0.0),
//                                cv.getAttributeAsDouble("MLEAF", 0.0),
//                                cv.getAttributeAsDouble("MLEAC", 0.0),
//                                cv.getAttributeAsDouble("FS", 0.0),
//                                cv.getAttributeAsInt("DP", 120),
//                                cv.getContig(),
//                                cv.getStart(),
//                                cv.getEnd(),
//                                cv.getReference().getDisplayString(),
//                                cv.getAlternateAllele(0).getDisplayString()
//                        ));
                    } else if (cv.isMNP()) {
                        if (cv.getReference().getBaseString().equals(SequenceUtils.reverseComplement(cv.getAlternateAllele(0).getBaseString()))) {
                            increment(stats, "INV", length, "FP");
                            usedV.add(cv);
                        } else {
                            increment(stats, "MNP", length, "FP");
                            usedV.add(cv);
                        }
                    } else if (cv.isSimpleInsertion() || (cv.getAlternateAllele(0).getDisplayString().contains(".") && cv.getAlternateAllele(0).length() > cv.getReference().length())) {
                        IndexedReference ref = cv.getContig().contains("HB3") ? REFERENCES.get("HB3") : REFERENCES.get("DD2");
                        String flankLeft = ref.getReferenceSequence().getSubsequenceAt(cv.getContig(), cv.getStart() - cv.getAlternateAllele(0).length(), cv.getEnd()).getBaseString();
                        String flankRight = ref.getReferenceSequence().getSubsequenceAt(cv.getContig(), cv.getStart(), cv.getEnd() + cv.getAlternateAllele(0).length()).getBaseString();

                        double sl = compareKmerSets(flankLeft, cv.getAlternateAllele(0).getDisplayString().replaceAll("\\.", ""), 11);
                        double sr = compareKmerSets(flankRight, cv.getAlternateAllele(0).getDisplayString().replaceAll("\\.", ""), 11);

                        if (sl > 0.5 || sr > 0.5) {
                            if (length > 20) {
                                increment(stats, "TD", length, "FP");
                                usedV.add(cv);
                            } else {
                                increment(stats, "STR_EXP", 0, "FP");
                                usedV.add(cv);
                            }
                        } else {
                            increment(stats, "INS", length, "FP");
                            usedV.add(cv);
                        }
                    } else if (cv.isSimpleDeletion() || (cv.getAlternateAllele(0).getDisplayString().contains(".") && cv.getAlternateAllele(0).length() < cv.getReference().length())) {
                        IndexedReference ref = cv.getContig().contains("HB3") ? REFERENCES.get("HB3") : REFERENCES.get("DD2");
                        String flankLeft = ref.getReferenceSequence().getSubsequenceAt(cv.getContig(), cv.getStart() - cv.getAlternateAllele(0).length(), cv.getEnd()).getBaseString();
                        String flankRight = ref.getReferenceSequence().getSubsequenceAt(cv.getContig(), cv.getStart(), cv.getEnd() + cv.getAlternateAllele(0).length()).getBaseString();

                        double sl = compareKmerSets(flankLeft, cv.getAlternateAllele(0).getDisplayString(), 11);
                        double sr = compareKmerSets(flankRight, cv.getAlternateAllele(0).getDisplayString(), 11);

                        if (length <= 20 && (sl > 0.5 || sr > 0.5 || flankLeft.contains(cv.getAlternateAllele(0).getDisplayString()) || flankRight.contains(cv.getAlternateAllele(0).getDisplayString()))) {
                            increment(stats, "STR_CON", 0, "FP");
                            usedV.add(cv);
                        } else {
                            increment(stats, "DEL", length, "FP");
                            usedV.add(cv);
                        }
                    } else if (cv.isIndel()) {
                        GeneratedVariant closestGv = closestSimulatedVariant(gvs, cv);
//                        log.info("indel");
//                    for (GeneratedVariant gv : gvs) {
//                        double score = computeConcreteSimilarity(gv, cv);
//                        if (score > 0.1) {
//                            score = computeConcreteSimilarity(gv, cv);
//                        }
//                    }
                    } else if (cv.getAlternateAllele(0).getDisplayString().contains(":") && !cv.getAlternateAllele(0).getDisplayString().contains("TANDEM")) {
                        if (cv.getAlternateAllele(0).getDisplayString().contains(cv.getContig())) {
                            increment(stats, "UNKNOWN", 0, "FP");
                        } else {
                            Interval l = convertAltToLocus(cv);
                            int j = 0;
                            increment(stats, "NAHR-INS", 0, "FP");
                        }
                    } else {
                        //int q = 0;
                        log.info("{}", cv);
                    }
                }

                if (!id1.equals(".")) { seenIds.add(id1); }
                if (!id2.equals(".")) { seenIds.add(id2); }
            }

            for (VariantContext cv : cvs) {
                if (!usedV.contains(cv)) {
                    unusedCVs.add(cv);
                }
            }
        }

        log.info("N: {} {} {}", allGvs.size(), allCvs.size(), unusedCVs.size());

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
//                String.format("%-7s", "C")
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
//                    tp + fn
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

        for (String key : Arrays.asList("UNKNOWN", "SNV", "STR_EXP", "STR_CON", "NAHR-INS", "NAHR-COMPLETE")) {
            Map<String, Integer> h = new LinkedHashMap<>();
            h.put("N", 0);
            h.put("TP", 0);
            h.put("FN", 0);
            h.put("FP", 0);

            stats.put(key, h);
        }

        for (String key : Arrays.asList("TD", "INV", "DEL", "MNP", "INS")) {
            for (String l : Arrays.asList("(1-100)", "(101-500)", "(501-1000)")) {
                Map<String, Integer> h = new LinkedHashMap<>();
                h.put("N", 0);
                h.put("TP", 0);
                h.put("FN", 0);
                h.put("FP", 0);

                stats.put(key + " " + l, h);
            }
        }

        return stats;
    }

    private Map<GeneratedVariant, Map<VariantContext, Double>> compareConcreteVariants(Set<GeneratedVariant> gvs, Set<VariantContext> cvs) {
        Map<GeneratedVariant, Map<VariantContext, Double>> ma = new HashMap<>();

        for (GeneratedVariant gv : gvs) {
            ma.put(gv, new HashMap<>());
            for (VariantContext cv : cvs) {
                if (!cv.isSymbolic() || cv.getAlternateAllele(0).getDisplayString().contains(".")) {
                    double score = computeConcreteSimilarity(gv, cv);

                    if (score > 0.5) {
                        ma.get(gv).put(cv, score);
                    }
                }
            }
        }

        return ma;
    }

    private Map<GeneratedVariant, Map<Pair<VariantContext, VariantContext>, Double>> compareSymbolicVariants(Set<GeneratedVariant> gvs, Set<VariantContext> cvs) {
        Map<GeneratedVariant, Map<Pair<VariantContext, VariantContext>, Double>> mb = new HashMap<>();

        for (GeneratedVariant gv : gvs) {
            mb.put(gv, new HashMap<>());
            for (VariantContext cv1 : cvs) {
                for (VariantContext cv2 : cvs) {
                    if (!cv1.equals(cv2)) {
                        double score = computeSymbolicSimilarity(gv, cv1, cv2);

                        if (score > 0.5) {
                            mb.get(gv).put(new Pair<>(cv1, cv2), computeSymbolicSimilarity(gv, cv1, cv2));
                        }
                    }
                }
            }
        }

        return mb;
    }

    private Interval convertAltToLocus(VariantContext cv) {
        String q = cv.getAlternateAllele(0).getDisplayString()
                .replaceAll("(HB3|DD2):", "")
                .replaceAll("[\\[\\]]", "")
                .replaceAll("^[ACGTN]+", "")
                .replaceAll("[ACGTN]+$", "")
                .replaceAll(":[+-]:.*", "");
                //.replaceAll("-.*", "");

        String[] qa = q.split("[:-]");
        String chr2 = qa[0];
        long start2 = Integer.parseInt(qa[1]);
        long end2 = start2 + 100;
        if (qa.length >= 3) {
            end2 = Integer.parseInt(qa[2]);
        }

        return new Interval(chr2, (int) start2, (int) end2);
    }

    private double computeConcreteSimilarity(GeneratedVariant gv, VariantContext cv) {
        String simulatedHaplotype = (gv.seedLeft + gv.newAllele + gv.seedRight).toUpperCase();
        String calledHaplotype = "";

        IndexedReference ref = cv.getContig().contains("HB3") ? REFERENCES.get("HB3") : REFERENCES.get("DD2");

        String hleft, hright;

        if (cv.getStart() == cv.getEnd()) {
            hleft = ref.getReferenceSequence().getSubsequenceAt(cv.getContig(), cv.getStart() - 100, cv.getStart() - 1).getBaseString();
            hright = ref.getReferenceSequence().getSubsequenceAt(cv.getContig(), cv.getEnd() + cv.getAlternateAllele(0).length() + 1, cv.getEnd() + cv.getAlternateAllele(0).length() + 100).getBaseString();
        } else {
            hleft = ref.getReferenceSequence().getSubsequenceAt(cv.getContig(), cv.getStart() - 100, cv.getStart() - 1).getBaseString();
            hright = ref.getReferenceSequence().getSubsequenceAt(cv.getContig(), cv.getEnd() + 1, cv.getEnd() + 100).getBaseString();
        }

        calledHaplotype = hleft + cv.getAlternateAllele(0).getDisplayString() + hright;

        return (!REQUIRE_CORRECT_PARENT || gv.start.getContig().equals(cv.getContig())) ? compareKmerSets(simulatedHaplotype, calledHaplotype, 11) : 0.0;
    }

    private double computeSymbolicSimilarity(GeneratedVariant gv, VariantContext cv1, VariantContext cv2) {
        double score = 0.0;
        if (gv.start.getContig().equals(cv1.getContig()) && gv.stop.getContig().equals(cv2.getContig())) {
            double leftDist = Math.min(Math.abs(gv.start.getStart() - cv1.getStart()), Math.abs(gv.start.getStart() - cv2.getStart()));
            double rightDist = Math.min(Math.abs(gv.stop.getStart() - cv2.getStart()), Math.abs(gv.stop.getStart() - cv1.getStart()));

            score = (5000.0 - (leftDist + rightDist)) / 5000.0;
        }

        return score;
    }

    private GeneratedVariant closestSimulatedVariant(Set<GeneratedVariant> gvs, VariantContext cv) {
        GeneratedVariant closestGv = null;
        int closestDist = Integer.MAX_VALUE;

        IndexedReference cvref = cv.getContig().contains("HB3") ? REFERENCES.get("HB3") : REFERENCES.get("DD2");
        int cvSeqIndex = cvref.getReferenceSequence().getSequenceDictionary().getSequenceIndex(cv.getContig());

        for (GeneratedVariant gv : gvs) {
            IndexedReference gvref1 = gv.start.getContig().contains("HB3") ? REFERENCES.get("HB3") : REFERENCES.get("DD2");
            int gvStartSeqIndex = gvref1.getReferenceSequence().getSequenceDictionary().getSequenceIndex(gv.start.getContig());

            IndexedReference gvref2 = gv.stop.getContig().contains("HB3") ? REFERENCES.get("HB3") : REFERENCES.get("DD2");
            int gvStopSeqIndex = gvref2.getReferenceSequence().getSequenceDictionary().getSequenceIndex(gv.start.getContig());

            if (cvSeqIndex == gvStartSeqIndex) {
                int dist = Math.abs(cv.getStart() - gv.start.getStart());

                if (dist < closestDist) {
                    closestGv = gv;
                    closestDist = dist;
                }
            }

            if (cvSeqIndex == gvStopSeqIndex) {
                int dist = Math.abs(cv.getStart() - gv.stop.getStart());

                if (dist < closestDist) {
                    closestGv = gv;
                    closestDist = dist;
                }
            }
        }

        return closestGv;
    }

    private VariantContext closestCalledVariant(GeneratedVariant gv, Set<VariantContext> cvs, boolean start) {
        VariantContext closestVc = null;
        int closestDist = Integer.MAX_VALUE;

        IndexedReference gvref1 = gv.start.getContig().contains("HB3") ? REFERENCES.get("HB3") : REFERENCES.get("DD2");
        int gvStartSeqIndex = gvref1.getReferenceSequence().getSequenceDictionary().getSequenceIndex(gv.start.getContig());

        IndexedReference gvref2 = gv.stop.getContig().contains("HB3") ? REFERENCES.get("HB3") : REFERENCES.get("DD2");
        int gvStopSeqIndex = gvref2.getReferenceSequence().getSequenceDictionary().getSequenceIndex(gv.start.getContig());

        for (VariantContext cv : cvs) {
            IndexedReference cvref = cv.getContig().contains("HB3") ? REFERENCES.get("HB3") : REFERENCES.get("DD2");
            int cvSeqIndex = cvref.getReferenceSequence().getSequenceDictionary().getSequenceIndex(cv.getContig());

            if (start) {
                boolean chrMatch = REQUIRE_CORRECT_PARENT ? cv.getContig().equalsIgnoreCase(gv.start.getContig()) : cvSeqIndex == gvStartSeqIndex;
                if (chrMatch) {
                    int dist = Math.abs(cv.getStart() - gv.start.getStart());

                    if (dist < closestDist && compatibleTypes(gv, cv)) {
                        closestVc = cv;
                        closestDist = dist;
                    }
                }
            } else {
                boolean chrMatch = REQUIRE_CORRECT_PARENT ? cv.getContig().equalsIgnoreCase(gv.stop.getContig()) : cvSeqIndex == gvStopSeqIndex;
                if (chrMatch) {
                    int dist = Math.abs(cv.getStart() - gv.stop.getStart());

                    if (dist < closestDist && compatibleTypes(gv, cv)) {
                        closestVc = cv;
                        closestDist = dist;
                    }
                }
            }
        }

        return closestDist < PROXIMAL_MATCH_RANGE ? closestVc : null;
    }

    private boolean compatibleTypes(GeneratedVariant gv, VariantContext vc) {
        if (gv.type.equals("SNV") && vc.isSNP()) { return true; }
        else if ((gv.type.equals("STR_EXP") || gv.type.equals("TD") || gv.type.equals("INS")) && (vc.isSimpleInsertion() || vc.isSymbolic())) { return true; }
        else if ((gv.type.equals("STR_CON") || gv.type.equals("DEL") || gv.type.equals("NAHR-INS")) && (vc.isSimpleDeletion() || vc.isSymbolic())) { return true; }
        else if ((gv.type.equals("INV") || gv.type.equals("MNP")) && (vc.isIndel() || vc.isMNP() || vc.isSymbolic())) { return true; }
        else if (gv.type.equals("NAHR-INS") && vc.isSymbolic()) { return true; }

        return false;
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
                    start = new Interval(refStartPieces[0], Integer.parseInt(refStartPieces[1]), Integer.parseInt(refStartPieces[2]), refStartPieces.length == 3, ".");

                    String[] refStopPieces = tr.get("refStop").split("[:-]");
                    stop = new Interval(refStopPieces[0], Integer.parseInt(refStopPieces[1]), Integer.parseInt(refStopPieces[2]), refStopPieces.length == 3, ".");
                } else {
                    start = new Interval(tr.get("refChr"), Integer.parseInt(tr.get("refStart")), Integer.parseInt(tr.get("refStart")));
                    stop = new Interval(tr.get("refChr"), Integer.parseInt(tr.get("refStop")), Integer.parseInt(tr.get("refStop")));
                }

                GeneratedVariant gv = new GeneratedVariant(
                        tr.get("type"),
                        Integer.parseInt(tr.get("chr")) - 1,
                        Integer.parseInt(tr.get("start")),
                        tr.get("old"),
                        tr.get("new"),
                        tr.get("sleft"),
                        tr.get("sright"),
                        tr.get("parent"),
                        start,
                        stop,
                        Integer.parseInt(tr.get("index"))
                );

                gvs.add(gv);
            }
        }

        return gvs;
    }

    private Set<VariantContext> loadCalledVariants(File variantFile) {
        VCFFileReader vcf = new VCFFileReader(variantFile, false);
        Set<VariantContext> vcs = new HashSet<>();

        IntervalTreeMap<VariantContext> it = new IntervalTreeMap<>();

        int nVariants = 0, nKeptVariants = 0;

        for (VariantContext vc : vcf) {
            nVariants++;

            IndexedReference ref = (vc.getContig().contains("HB3")) ? REFERENCES.get("HB3") : REFERENCES.get("DD2");
            IndexedReference aref = (vc.getContig().contains("HB3")) ? REFERENCES.get("DD2") : REFERENCES.get("HB3");

            if (!vc.getContig().contains("contig") &&
                !vc.getAlternateAllele(0).getDisplayString().contains("contig") &&
                vc.getEnd() + vc.getAlternateAllele(0).length() + 100 < ref.getReferenceSequence().getSequence(vc.getContig()).length() &&
                vc.getEnd() - vc.getStart() < 10000 &&
                vc.getAttributeAsInt("flankMappingQuality", 60) > 0 &&
                vc.getAttributeAsInt("DP", 120) > 80 &&
                vc.getAttributeAsInt("DP", 120) < 200 &&
                vc.getAttributeAsDouble("SOR", 0.0) < 3.0 &&
                vc.getAttributeAsDouble("FS", 0.0) < 1.0 &&
                vc.getAttributeAsDouble("MQ", 60.0) > 9.0 &&
                (!EXCLUDE_DUPLICATES_IN_SYNTENIC_LOCI || !it.containsOverlapping(new Interval(vc.getContig(), vc.getStart(), vc.getEnd())))
            ) {
                VariantContext vcnew = vc;

                if (vc.isSymbolic()) {
                    if (vc.hasAttribute("CONSENSUS")) {
                        vcnew = new VariantContextBuilder(vc)
                            .alleles(Arrays.asList(vc.getReference().getBaseString(), vc.getAttributeAsString("CONSENSUS", "N")))
                            .make();
                    } else {
                        String refAllele = ref.getReferenceSequence().getSubsequenceAt(vc.getContig(), vc.getStart(), vc.getEnd()).getBaseString();
                        String altAllele = null;
                        if (vc.getAlternateAllele(0).getDisplayString().contains("DEL")) {
                            altAllele = refAllele.substring(0, 1);
                        } else if (vc.getAlternateAllele(0).getDisplayString().contains("INV")) {
                            altAllele = SequenceUtils.reverseComplement(refAllele.substring(1));
                        } else if (vc.getAlternateAllele(0).getDisplayString().contains("DUP")) {
                            altAllele = refAllele + refAllele;
                        }

                        if (altAllele != null) {
                            vcnew = new VariantContextBuilder(vc)
                                    .alleles(Arrays.asList(refAllele, altAllele))
                                    .make();
                        }
                    }
                }

                List<SAMRecord> srs = new ArrayList<>();
                if (!vcnew.isSymbolic() && vcnew.getEnd() + vcnew.getAlternateAllele(0).length() + 100 < ref.getReferenceSequence().getSequence(vcnew.getContig()).length()) {
                    String left = ref.getReferenceSequence().getSubsequenceAt(vcnew.getContig(), vcnew.getStart() - 50, vcnew.getStart() - 1).getBaseString();
                    String right = ref.getReferenceSequence().getSubsequenceAt(vcnew.getContig(), vcnew.getEnd() + vcnew.getAlternateAllele(0).length() + 1, vcnew.getEnd() + vcnew.getAlternateAllele(0).length() + 50).getBaseString();
                    String hap = left + vcnew.getAlternateAllele(0).getDisplayString() + right;

                    srs = aref.align(hap);
                    srs.removeIf(s -> (s.getMappingQuality() < 10 || !s.getCigarString().equals(hap.length() + "M")));
                }

                if (srs.size() == 0) {
                    vcs.add(vcnew);
                    nKeptVariants++;
                }
            }

            if (!vc.getContig().contains("contig") && vc.getStart() + 100 < ref.getReferenceSequence().getSequence(vc.getContig()).length()) {
                String refHaplotype = ref.getReferenceSequence().getSubsequenceAt(vc.getContig(), vc.getStart() - 100, vc.getStart() + 100).getBaseString();
                List<SAMRecord> altAlignments = aref.align(refHaplotype);
                altAlignments.removeIf(s -> (s.getMappingQuality() < 10));
                if (altAlignments.size() == 1) {
                    int refIndex = ref.getReferenceSequence().getSequenceDictionary().getSequenceIndex(vc.getContig());
                    int altIndex = aref.getReferenceSequence().getSequenceDictionary().getSequenceIndex(altAlignments.get(0).getContig());

                    if (refIndex == altIndex) {
                        it.put(new Interval(altAlignments.get(0).getContig(), altAlignments.get(0).getAlignmentStart(), altAlignments.get(0).getAlignmentEnd()), vc);
                    }
                }
            }
        }

        log.debug("Loaded {} variants, retained {} variants", nVariants, nKeptVariants);

        return vcs;
    }
}

