package uk.ac.ox.well.cortexjdk.playground.eval;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.commands.simulate.generators.GeneratedVariant;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.containers.ContainerUtils;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class EvaluateVariants extends Module {
    @Argument(fullName="truth", shortName="t", doc="Truth table")
    public ArrayList<File> TRUTH_TABLES;

    @Argument(fullName="novels", shortName="n", doc="Novels")
    public ArrayList<File> NOVELS;

    @Argument(fullName="calls", shortName="c", doc="Calls")
    public ArrayList<File> CALLS;

    @Argument(fullName="references", shortName="R", doc="Reference(s)", required=false)
    public HashMap<String, IndexedReference> REFERENCES;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Output
    public PrintStream out;

    private class Validation {
        public boolean correctChromosome = false;
        public boolean correctStart = false;
        public boolean correctStop = false;
        public boolean correctType = false;
        public boolean correctRefAllele = false;
        public boolean correctAltAllele = false;
        public boolean correctHaplotype = false;
        //public boolean correctPhase = false;
        public boolean closeToStart = false;
        public boolean closeToStop = false;
        public boolean breakpointWithinRange = false;
        public boolean correctLength = false;
        public boolean isSymbolic = false;
        public Set<VariantContext> vcs = new HashSet<>();

        @Override
        public String toString() {
            return "Validation{" +
                    "correctChromosome=" + correctChromosome +
                    ", correctStart=" + correctStart +
                    ", correctStop=" + correctStop +
                    ", correctType=" + correctType +
                    ", correctRefAllele=" + correctRefAllele +
                    ", correctAltAllele=" + correctAltAllele +
                    ", correctHaplotype=" + correctHaplotype +
                    //", correctPhase=" + correctPhase +
                    ", closeToStart=" + closeToStart +
                    ", closeToStop=" + closeToStop +
                    ", breakpointWithinRange=" + breakpointWithinRange +
                    ", isSymbolic=" + isSymbolic +
                    ", vcs=" + vcs +
                    '}';
        }
    }

    @Override
    public void execute() {
        evaluateCalls(TRUTH_TABLES, CALLS, NOVELS);
    }

    private void evaluateCalls(List<File> truthTableF, List<File> callsF, List<File> novelsF) {
        Map<String, Map<String, Integer>> simTable = new LinkedHashMap<>();
        for (String id : Arrays.asList(
                "SNVs",
                "Insertions (random, 1-100)",
                "Insertions (random, 101-500)",
                "Insertions (random, 501-1000)",
                "Insertions (tandem duplications, 1-100)",
                "Insertions (tandem duplications, 101-500)",
                "Insertions (tandem duplications, 501-1000)",
                "Insertions (STR expansions)",
                "Deletions (random, 1-100)",
                "Deletions (random, 101-500)",
                "Deletions (random, 501-1000)",
                "Deletions (STR contractions)",
                "MNVs (random, 1-100)",
                "MNVs (random, 101-500)",
                "MNVs (random, 501-1000)",
                "MNVs (inversions, 1-100)",
                "MNVs (inversions, 101-500)",
                "MNVs (inversions, 501-1000)",
                "NAHRs",
                "NAHRs (complete)"
        )) {
            simTable.put(id, new LinkedHashMap<>());
            simTable.get(id).put("N", 0);
            simTable.get(id).put("Novels", 0);
            simTable.get(id).put("Novels recovered", 0);
            simTable.get(id).put("TP", 0);
            simTable.get(id).put("FP", 0);
            simTable.get(id).put("FN", 0);
        }

        int numSimulatedVariants = 0;
        int numSimulatedVariants2 = 0;

        for (int i = 0; i < truthTableF.size(); i++) {
            File truthTable = truthTableF.get(i);
            File callsFile = callsF.get(i);
            File novelsFile = novelsF.get(i);

            if (!truthTable.getAbsolutePath().split("/")[8].equals(callsFile.getAbsolutePath().split("/")[8]) ||
                !truthTable.getAbsolutePath().split("/")[8].equals(novelsFile.getAbsolutePath().split("/")[8])) {
                throw new CortexJDKException("Mismatch sample names: " +
                        truthTable.getAbsolutePath().split("/")[8] + " " +
                        callsFile.getAbsolutePath().split("/")[8] + " " +
                        novelsFile.getAbsolutePath().split("/")[8]);
            }

            TableReader tna = new TableReader(novelsFile, "index", "len", "i", "kmer");
            Map<Integer, Integer> indexToNumNovels = new HashMap<>();
            for (Map<String, String> te : tna) {
                int index = Integer.valueOf(te.get("index"));
                ContainerUtils.increment(indexToNumNovels, index);
            }

            List<GeneratedVariant> vs = new ArrayList<>();

            TableReader truthTr = new TableReader(truthTable);

            Map<Integer, List<Map<String, String>>> recombs = new TreeMap<>();
            for (Map<String, String> tr : truthTr) {
                if (tr.get("type").equals("RECOMB")) {
                    int chrIndex = Integer.valueOf(tr.get("chr"));
                    if (!recombs.containsKey(chrIndex)) {
                        recombs.put(chrIndex, new ArrayList<>());
                    }

                    recombs.get(chrIndex).add(tr);
                } else {
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

                    vs.add(gv);
                    numSimulatedVariants++;
                }
            }

            //List<Pair<List<String>, List<Integer>>> seqs = recombine(recombs);

            VCFFileReader vcf = new VCFFileReader(callsFile, false);
            Set<VariantContext> vcs = new HashSet<>();
            for (VariantContext vc : vcf) {
                vcs.add(vc);
            }

            //List<GeneratedVariant> gvs = updateLociOfSimulatedVariants(vs);
            Map<GeneratedVariant, Set<VariantContext>> cvs = clusterCalledVariants(vcs, vs);
            //Map<GeneratedVariant, Set<VariantContext>> rvs = filterNAHRCalls(vcs, gvs, cvs);

            numSimulatedVariants2 += cvs.size() - 1;

            for (GeneratedVariant gv : cvs.keySet()) {
                if (!gv.type.equals("rejected")) {
                    Set<VariantContext> toRemove = new HashSet<>();

                    for (VariantContext v : cvs.get(gv)) {
                        Validation val = new Validation();

                        //String bg = v.getAttributeAsString("BACKGROUND", REFERENCES.keySet().iterator().next());
                        String bg = v.getContig().contains("HB3") ? "HB3" : "DD2";
                        int seqIndex = REFERENCES.get(bg).getReferenceSequence().getSequenceDictionary().getSequenceIndex(v.getContig());

                        if (gv.seqIndex + 1 == seqIndex) { val.correctChromosome = true; }
                        if (gv.getPosIndex() + 1 == v.getStart()) { val.correctStart = true; }
                        if (gv.getPosIndex() + 1 == v.getEnd()) { val.correctStop = true; }
                        if (gv.getOldAllele().equalsIgnoreCase(v.getReference().getBaseString())) { val.correctRefAllele = true; }
                        if (gv.getNewAllele().equalsIgnoreCase(v.getAltAlleleWithHighestAlleleCount().getBaseString())) { val.correctAltAllele = true; }
                        if (isSameType(gv, v)) { val.correctType = true; }
                        //if (gv.parent.equalsIgnoreCase(v.getAttributeAsString("BACKGROUND", ""))) { val.correctPhase = true; }
                        if (Math.abs(gv.start.getStart() - v.getStart()) < 100) { val.closeToStart = true; }
                        if (Math.abs(gv.stop.getStart() - v.getStart() - v.getAltAlleleWithHighestAlleleCount().length()) < 100) { val.closeToStop = true; }
                        if (!gv.start.getContig().equals(gv.stop.getContig())) {
                            String[] p = v.getAltAlleleWithHighestAlleleCount().getDisplayString().split(":");
                            if (p.length == 5) {
                                String[] l = p[2].split("-");
                                Interval ql = new Interval(p[1], Integer.valueOf(l[0]), Integer.valueOf(l[1]), p[3].equals("-"), ".");
                                if (ql.intersects(gv.stop)) {
                                    val.breakpointWithinRange = true;
                                }
                            }
                        }
//                        if ((gv.getType().contains("NAHR")) || (gv.start.getContig().equals(gv.stop.getContig()) && Math.abs(gv.stop.getStart() - gv.start.getStart()) > 100)) {
//                            if (v.isSymbolic()) {
//                                val.isSymbolic = true;
//                            }
//                        }

                        String generatedChildHap = (gv.seedLeft + gv.newAllele + gv.seedRight).toUpperCase();
                        String calledChildHap = (v.getAttributeAsString("CHILD_HAP", "")).toUpperCase();

                        if (generatedChildHap.contains(calledChildHap)) { val.correctHaplotype = true; }

                        float expLength = (float) (gv.newAllele.length() - gv.getOldAllele().length());
                        float actLength = (float) (v.getAltAlleleWithHighestAlleleCount().length() - v.getReference().length());
                        if (0.80*actLength <= expLength && expLength <= 1.20*actLength) {
                            val.correctLength = true;
                        }

                        if (gv.type.equals("SNV")) {
                            if (val.correctChromosome && val.correctType && val.correctRefAllele && val.correctAltAllele && ((val.correctStart && val.correctStop) || (val.closeToStart && val.closeToStop))) {
                            } else {
                                toRemove.add(v);
                            }
                        } else if (gv.type.contains("NAHR-INS")) {
                            int score = 0;
                            score += val.correctType ? 1 : 0;
                            score += val.correctStart ? 1 : 0;
                            score += val.correctStop ? 1 : 0;
                            score += val.closeToStop ? 1 : 0;
                            score += val.closeToStart ? 1 : 0;
                            //score += val.correctPhase ? 1 : 0;
                            score += val.breakpointWithinRange ? 1 : 0;
                            if (val.correctChromosome && (val.correctHaplotype || score > 2)) {
                                //log.info("");
                            } else {
                                toRemove.add(v);
                            }
                        } else {
                            if (v.isSymbolic()) {
                                int score = 0;
                                score += val.correctStart ? 1 : 0;
                                score += val.correctStop ? 1 : 0;
                                score += val.closeToStop ? 1 : 0;
                                score += val.closeToStart ? 1 : 0;
                                //score += val.correctPhase ? 1 : 0;

                                if (val.correctChromosome && (val.correctHaplotype || score >= 2)) {
                                } else {
                                    toRemove.add(v);
                                }
                            } else {
                                int score = 0;
                                score += val.correctType ? 1 : 0;
                                score += val.correctLength ? 1 : 0;
                                score += val.correctStart ? 1 : 0;
                                score += val.correctStop ? 1 : 0;
                                score += val.closeToStop ? 1 : 0;
                                score += val.closeToStart ? 1 : 0;

                                if (val.correctChromosome && (val.correctHaplotype || score >= 2)) {
                                } else {
                                    toRemove.add(v);
                                }
                            }
                        }
                    }

                    cvs.get(gv).removeAll(toRemove);

                    //log.info("");
                }
            }

            Set<VariantContext> unused = new HashSet<>();
            for (VariantContext v : vcs) {
                boolean found = false;
                for (GeneratedVariant gv : cvs.keySet()) {
                    if (!gv.type.equals("rejected")) {
                        for (VariantContext v2 : cvs.get(gv)) {
                            if (v.getAttributeAsString("CALL_ID", "0").equals(v2.getAttributeAsString("CALL_ID", "0"))) {
                                found = true;
                                break;
                            }
                        }
                    } else {
                        //log.info("");
                    }
                }

                if (!found && !v.isSymbolic()) {
                    unused.add(v);
                }
            }

            /*
            "SNVs",
            "Insertions (random, 1-100)",
            "Insertions (random, 101-500)",
            "Insertions (random, 501-1000)",
            "Insertions (STR expansions)",
            "Insertions (tandem duplications)",
            "Deletions (random, 1-100)",
            "Deletions (random, 101-500)",
            "Deletions (random, 501-1000)",
            "Deletions (STR contractions)",
            "MNVs (random, 1-100)",
            "MNVs (random, 101-500)",
            "MNVs (random, 501-1000)",
            "MNVs (inversions, 1-100)",
            "MNVs (inversions, 101-500)",
            "MNVs (inversions, 501-1000)",
            "NAHRs"
            */

            Map<String, Set<VariantContext>> vcsByPartition = new HashMap<>();
            for (VariantContext vc : vcs) {
                if (!vcsByPartition.containsKey(vc.getAttributeAsString("PARTITION_NAME", ""))) {
                    vcsByPartition.put(vc.getAttributeAsString("PARTITION_NAME", ""), new HashSet<>());
                }

                vcsByPartition.get(vc.getAttributeAsString("PARTITION_NAME", "")).add(vc);
            }

            Set<String> partitionFilter = new HashSet<>();
            for (String partitionName : vcsByPartition.keySet()) {
                Set<String> chroms = new HashSet<>();
                for (VariantContext v : vcsByPartition.get(partitionName)) {
                    chroms.add(v.getContig());
                    if (v.getAltAlleleWithHighestAlleleCount().getDisplayString().contains(":")) {
                        chroms.add(v.getAltAlleleWithHighestAlleleCount().getDisplayString().split(":")[1]);
                    }
                }

                if (chroms.size() > 2) {
                    partitionFilter.add(partitionName);
                }
            }

            int nahrs_exp = 0;
            int nahrs_act = 0;
            for (GeneratedVariant gv : cvs.keySet()) {
                if (gv.type.equals("rejected")) {
                    for (VariantContext vc : cvs.get(gv)) {
                        if (vc.getContig().contains("contig") || vc.getContig().contains("partition")) {
                            partitionFilter.add(vc.getAttributeAsString("PARTITION_NAME", "unknown"));
                        }
                    }

                    for (VariantContext vc : cvs.get(gv)) {
                        if (partitionFilter.contains(vc.getAttributeAsString("PARTITION_NAME", "unknown"))) {
                            // filtered out
                        } else if (vc.isSNP()) {
                            ContainerUtils.increment(simTable.get("SNVs"), "FP");
                        } else if (vc.isSimpleInsertion()) {
                            if (vc.getAltAlleleWithHighestAlleleCount().length() <= 100) {
                                ContainerUtils.increment(simTable.get("Insertions (random, 1-100)"), "FP");
                            } else if (vc.getAltAlleleWithHighestAlleleCount().length() <= 500) {
                                ContainerUtils.increment(simTable.get("Insertions (random, 101-500)"), "FP");
                            } else {
                                ContainerUtils.increment(simTable.get("Insertions (random, 501-1000)"), "FP");
                            }
                        } else if (vc.isSimpleDeletion()) {
                            if (vc.getReference().length() <= 100) {
                                ContainerUtils.increment(simTable.get("Deletions (random, 1-100)"), "FP");
                            } else if (vc.getReference().length() <= 500) {
                                ContainerUtils.increment(simTable.get("Deletions (random, 101-500)"), "FP");
                            } else {
                                ContainerUtils.increment(simTable.get("Deletions (random, 501-1000)"), "FP");
                            }
                        } else if (vc.isMNP()) {
                            if (vc.getReference().length() <= 100) {
                                ContainerUtils.increment(simTable.get("MNVs (random, 1-100)"), "FP");
                            } else if (vc.getReference().length() <= 500) {
                                ContainerUtils.increment(simTable.get("MNVs (random, 101-500)"), "FP");
                            } else {
                                ContainerUtils.increment(simTable.get("MNVs (random, 501-1000)"), "FP");
                            }
                        } else if (vc.isSymbolic() && vc.getAltAlleleWithHighestAlleleCount().getDisplayString().contains(":")) {
                            if (!vc.getContig().equals(vc.getAltAlleleWithHighestAlleleCount().getDisplayString().split(":")[1])) {
                                ContainerUtils.increment(simTable.get("NAHRs"), "FP");
                            }
                        }
                    }
                } else if (gv.type.equals("SNV")) {
                    ContainerUtils.increment(simTable.get("SNVs"), "N");
                    if (indexToNumNovels.containsKey(gv.index)) simTable.get("SNVs").put("Novels", simTable.get("SNVs").get("Novels") + indexToNumNovels.get(gv.index));
                    if (cvs.get(gv).size() > 0) {
                        ContainerUtils.increment(simTable.get("SNVs"), "TP");
                        if (indexToNumNovels.containsKey(gv.index)) simTable.get("SNVs").put("Novels recovered", simTable.get("SNVs").get("Novels recovered") + indexToNumNovels.get(gv.index));
                    }
                    else {
                        ContainerUtils.increment(simTable.get("SNVs"), "FN");
                    }
                } else if (gv.type.equals("INS")) {
                    if (gv.getNewAllele().length() <= 100) {
                        ContainerUtils.increment(simTable.get("Insertions (random, 1-100)"), "N");
                        if (indexToNumNovels.containsKey(gv.index)) simTable.get("Insertions (random, 1-100)").put("Novels", simTable.get("Insertions (random, 1-100)").get("Novels") + indexToNumNovels.get(gv.index));
                        if (cvs.get(gv).size() > 0) {
                            ContainerUtils.increment(simTable.get("Insertions (random, 1-100)"), "TP");
                            if (indexToNumNovels.containsKey(gv.index)) simTable.get("Insertions (random, 1-100)").put("Novels recovered", simTable.get("Insertions (random, 1-100)").get("Novels recovered") + indexToNumNovels.get(gv.index));
                        }
                        else { ContainerUtils.increment(simTable.get("Insertions (random, 1-100)"), "FN"); }
                    } else if (gv.getNewAllele().length() <= 500) {
                        ContainerUtils.increment(simTable.get("Insertions (random, 101-500)"), "N");
                        if (indexToNumNovels.containsKey(gv.index)) simTable.get("Insertions (random, 101-500)").put("Novels", simTable.get("Insertions (random, 101-500)").get("Novels") + indexToNumNovels.get(gv.index));
                        if (cvs.get(gv).size() > 0) {
                            ContainerUtils.increment(simTable.get("Insertions (random, 101-500)"), "TP");
                            if (indexToNumNovels.containsKey(gv.index)) simTable.get("Insertions (random, 101-500)").put("Novels recovered", simTable.get("Insertions (random, 101-500)").get("Novels recovered") + indexToNumNovels.get(gv.index));
                        }
                        else { ContainerUtils.increment(simTable.get("Insertions (random, 101-500)"), "FN"); }
                    } else {
                        ContainerUtils.increment(simTable.get("Insertions (random, 501-1000)"), "N");
                        if (indexToNumNovels.containsKey(gv.index)) simTable.get("Insertions (random, 501-1000)").put("Novels", simTable.get("Insertions (random, 501-1000)").get("Novels") + indexToNumNovels.get(gv.index));
                        if (cvs.get(gv).size() > 0) {
                            ContainerUtils.increment(simTable.get("Insertions (random, 501-1000)"), "TP");
                        }
                        else { ContainerUtils.increment(simTable.get("Insertions (random, 501-1000)"), "FN"); }
                    }
                } else if (gv.type.equals("STR_EXP")) {
                    ContainerUtils.increment(simTable.get("Insertions (STR expansions)"), "N");
                    if (indexToNumNovels.containsKey(gv.index)) simTable.get("Insertions (STR expansions)").put("Novels", simTable.get("Insertions (STR expansions)").get("Novels") + indexToNumNovels.get(gv.index));
                    if (cvs.get(gv).size() > 0) {
                        ContainerUtils.increment(simTable.get("Insertions (STR expansions)"), "TP");
                    }
                    else { ContainerUtils.increment(simTable.get("Insertions (STR expansions)"), "FN"); }
                } else if (gv.type.equals("TD")) {
                    if (Math.abs(gv.getNewAllele().length() - gv.getOldAllele().length()) <= 100) {
                        ContainerUtils.increment(simTable.get("Insertions (tandem duplications, 1-100)"), "N");
                        if (indexToNumNovels.containsKey(gv.index)) simTable.get("Insertions (tandem duplications, 1-100)").put("Novels", simTable.get("Insertions (tandem duplications, 1-100)").get("Novels") + indexToNumNovels.get(gv.index));
                        if (cvs.get(gv).size() > 0) {
                            ContainerUtils.increment(simTable.get("Insertions (tandem duplications, 1-100)"), "TP");
                        } else { ContainerUtils.increment(simTable.get("Insertions (tandem duplications, 1-100)"), "FN"); }
                    } else if (Math.abs(gv.getNewAllele().length() - gv.getOldAllele().length()) <= 500) {
                        ContainerUtils.increment(simTable.get("Insertions (tandem duplications, 101-500)"), "N");
                        if (indexToNumNovels.containsKey(gv.index)) simTable.get("Insertions (tandem duplications, 101-500)").put("Novels", simTable.get("Insertions (tandem duplications, 101-500)").get("Novels") + indexToNumNovels.get(gv.index));
                        if (cvs.get(gv).size() > 0) {
                            ContainerUtils.increment(simTable.get("Insertions (tandem duplications, 101-500)"), "TP");
                        } else { ContainerUtils.increment(simTable.get("Insertions (tandem duplications, 101-500)"), "FN"); }
                    } else {
                        ContainerUtils.increment(simTable.get("Insertions (tandem duplications, 501-1000)"), "N");
                        if (indexToNumNovels.containsKey(gv.index)) simTable.get("Insertions (tandem duplications, 501-1000)").put("Novels", simTable.get("Insertions (tandem duplications, 501-1000)").get("Novels") + indexToNumNovels.get(gv.index));
                        if (cvs.get(gv).size() > 0) {
                            ContainerUtils.increment(simTable.get("Insertions (tandem duplications, 501-1000)"), "TP");
                        } else { ContainerUtils.increment(simTable.get("Insertions (tandem duplications, 501-1000)"), "FN"); }
                    }
                } else if (gv.type.equals("DEL") || gv.type.equals("NAHR-DEL")) {
                    if (gv.type.equals("NAHR-DEL")) {
                        ContainerUtils.increment(simTable.get("NAHRs (complete)"), "N");
                    }

                    if (gv.getOldAllele().length() <= 100) {
                        ContainerUtils.increment(simTable.get("Deletions (random, 1-100)"), "N");
                        if (indexToNumNovels.containsKey(gv.index)) simTable.get("Deletions (random, 1-100)").put("Novels", simTable.get("Deletions (random, 1-100)").get("Novels") + indexToNumNovels.get(gv.index));
                        if (cvs.get(gv).size() > 0) {
                            ContainerUtils.increment(simTable.get("Deletions (random, 1-100)"), "TP");
                        }
                        else { ContainerUtils.increment(simTable.get("Deletions (random, 1-100)"), "FN"); }
                    } else if (gv.getOldAllele().length() <= 500) {
                        ContainerUtils.increment(simTable.get("Deletions (random, 101-500)"), "N");
                        if (indexToNumNovels.containsKey(gv.index)) simTable.get("Deletions (random, 101-500)").put("Novels", simTable.get("Deletions (random, 101-500)").get("Novels") + indexToNumNovels.get(gv.index));
                        if (cvs.get(gv).size() > 0) {
                            ContainerUtils.increment(simTable.get("Deletions (random, 101-500)"), "TP");
                        }
                        else { ContainerUtils.increment(simTable.get("Deletions (random, 101-500)"), "FN"); }
                    } else {
                        ContainerUtils.increment(simTable.get("Deletions (random, 501-1000)"), "N");
                        if (indexToNumNovels.containsKey(gv.index)) simTable.get("Deletions (random, 501-1000)").put("Novels", simTable.get("Deletions (random, 501-1000)").get("Novels") + indexToNumNovels.get(gv.index));
                        if (cvs.get(gv).size() > 0) {
                            ContainerUtils.increment(simTable.get("Deletions (random, 501-1000)"), "TP");
                        }
                        else { ContainerUtils.increment(simTable.get("Deletions (random, 501-1000)"), "FN"); }
                    }
                } else if (gv.type.equals("STR_CON")) {
                    ContainerUtils.increment(simTable.get("Deletions (STR contractions)"), "N");
                    if (indexToNumNovels.containsKey(gv.index)) simTable.get("Deletions (STR contractions)").put("Novels", simTable.get("Deletions (STR contractions)").get("Novels") + indexToNumNovels.get(gv.index));
                    if (cvs.get(gv).size() > 0) {
                        ContainerUtils.increment(simTable.get("Deletions (STR contractions)"), "TP");
                    }
                    else { ContainerUtils.increment(simTable.get("Deletions (STR contractions)"), "FN"); }
                } else if (gv.type.equals("MNP")) {
                    if (gv.getOldAllele().length() <= 100) {
                        ContainerUtils.increment(simTable.get("MNVs (random, 1-100)"), "N");
                        if (indexToNumNovels.containsKey(gv.index)) simTable.get("MNVs (random, 1-100)").put("Novels", simTable.get("MNVs (random, 1-100)").get("Novels") + indexToNumNovels.get(gv.index));
                        if (cvs.get(gv).size() > 0) {
                            ContainerUtils.increment(simTable.get("MNVs (random, 1-100)"), "TP");
                        }
                        else { ContainerUtils.increment(simTable.get("MNVs (random, 1-100)"), "FN"); }
                    } else if (gv.getOldAllele().length() <= 500) {
                        ContainerUtils.increment(simTable.get("MNVs (random, 101-500)"), "N");
                        if (indexToNumNovels.containsKey(gv.index)) simTable.get("MNVs (random, 101-500)").put("Novels", simTable.get("MNVs (random, 101-500)").get("Novels") + indexToNumNovels.get(gv.index));
                        if (cvs.get(gv).size() > 0) {
                            ContainerUtils.increment(simTable.get("MNVs (random, 101-500)"), "TP");
                        }
                        else { ContainerUtils.increment(simTable.get("MNVs (random, 101-500)"), "FN"); }
                    } else {
                        ContainerUtils.increment(simTable.get("MNVs (random, 501-1000)"), "N");
                        if (indexToNumNovels.containsKey(gv.index)) simTable.get("MNVs (random, 501-1000)").put("Novels", simTable.get("MNVs (random, 501-1000)").get("Novels") + indexToNumNovels.get(gv.index));
                        if (cvs.get(gv).size() > 0) {
                            ContainerUtils.increment(simTable.get("MNVs (random, 501-1000)"), "TP");
                        }
                        else { ContainerUtils.increment(simTable.get("MNVs (random, 501-1000)"), "FN"); }
                    }
                } else if (gv.type.equals("INV")) {
                    if (gv.getOldAllele().length() <= 100) {
                        ContainerUtils.increment(simTable.get("MNVs (inversions, 1-100)"), "N");
                        if (indexToNumNovels.containsKey(gv.index)) simTable.get("MNVs (inversions, 1-100)").put("Novels", simTable.get("MNVs (inversions, 1-100)").get("Novels") + indexToNumNovels.get(gv.index));
                        if (cvs.get(gv).size() > 0) {
                            ContainerUtils.increment(simTable.get("MNVs (inversions, 1-100)"), "TP");
                        }
                        else { ContainerUtils.increment(simTable.get("MNVs (inversions, 1-100)"), "FN"); }
                    } else if (gv.getOldAllele().length() <= 500) {
                        ContainerUtils.increment(simTable.get("MNVs (inversions, 101-500)"), "N");
                        if (indexToNumNovels.containsKey(gv.index)) simTable.get("MNVs (inversions, 101-500)").put("Novels", simTable.get("MNVs (inversions, 101-500)").get("Novels") + indexToNumNovels.get(gv.index));
                        if (cvs.get(gv).size() > 0) {
                            ContainerUtils.increment(simTable.get("MNVs (inversions, 101-500)"), "TP");
                        }
                        else { ContainerUtils.increment(simTable.get("MNVs (inversions, 101-500)"), "FN"); }
                    } else {
                        ContainerUtils.increment(simTable.get("MNVs (inversions, 501-1000)"), "N");
                        if (indexToNumNovels.containsKey(gv.index)) simTable.get("MNVs (inversions, 501-1000)").put("Novels", simTable.get("MNVs (inversions, 501-1000)").get("Novels") + indexToNumNovels.get(gv.index));
                        if (cvs.get(gv).size() > 0) {
                            ContainerUtils.increment(simTable.get("MNVs (inversions, 501-1000)"), "TP");
                        }
                        else { ContainerUtils.increment(simTable.get("MNVs (inversions, 501-1000)"), "FN"); }
                    }
                } else if (gv.type.equals("NAHR-INS")) {
                    nahrs_exp++;
                    ContainerUtils.increment(simTable.get("NAHRs"), "N");
                    if (indexToNumNovels.containsKey(gv.index)) simTable.get("NAHRs").put("Novels", simTable.get("NAHRs").get("Novels") + indexToNumNovels.get(gv.index));
                    if (indexToNumNovels.containsKey(gv.index)) simTable.get("NAHRs (complete)").put("Novels", simTable.get("NAHRs (complete)").get("Novels") + indexToNumNovels.get(gv.index));

                    if (cvs.get(gv).size() >= 1) {
                        ContainerUtils.increment(simTable.get("NAHRs"), "TP");
                        nahrs_act++;
                    } else {
                        ContainerUtils.increment(simTable.get("NAHRs"), "FN");
                    }
                } else {
                    log.info("??? {} gv");
                }
            }

//            if (nahrs_exp > 0) {
//                ContainerUtils.increment(simTable.get("NAHRs (complete)"), "N");
//
//                if (nahrs_act == nahrs_exp) {
//                    ContainerUtils.increment(simTable.get("NAHRs (complete)"), "TP");
//                } else {
//                    ContainerUtils.increment(simTable.get("NAHRs (complete)"), "FN");
//                }
//            }

            //log.info("NAHR: {} {}", nahrs_exp, nahrs_act);
        }

        List<String> columns = new ArrayList<>();
        columns.add(String.format("%-43s", "type"));
        for (String column : simTable.get(simTable.keySet().iterator().next()).keySet()) {
            String format = "%-" + (column.length() + 2) + "s";
            columns.add(String.format(format, column));
        }
        columns.add("F1");
        out.println(Joiner.on("\t").join(columns));

        for (String type : simTable.keySet()) {
            List<String> row = new ArrayList<>();
            for (String column : simTable.get(type).keySet()) {
                String format = "%-" + (column.length() + 2) + "s";
                row.add(String.format(format, simTable.get(type).get(column)));
            }

            double tp = (double) simTable.get(type).get("TP");
            double fp = (double) simTable.get(type).get("FP");
            double fn = (double) simTable.get(type).get("FN");
            double pr = (tp == 0 && fp == 0 && fn == 0) ? 0.0 : tp / (tp + fp);
            double re = (tp == 0 && fp == 0 && fn == 0) ? 0.0 : tp / (tp + fn);
            double f1 = (tp == 0 && fp == 0 && fn == 0) ? 0.0 : (2*tp) / ((2*tp) + fp + fn);

            out.println(String.format("%-43s", type) + "\t" + Joiner.on("\t").join(row) + "\t" + String.format("%.2f", f1));
        }
    }

    private Map<GeneratedVariant, Set<VariantContext>> filterNAHRCalls(Set<VariantContext> vcs, List<GeneratedVariant> gvs, Map<GeneratedVariant, Set<VariantContext>> clusteredVariants) {
        Map<GeneratedVariant, Set<VariantContext>> filteredVariants = new HashMap<>();

        GeneratedVariant reject = new GeneratedVariant("rejected", 0, 0, "", "");

        for (GeneratedVariant gv : clusteredVariants.keySet()) {
            Set<VariantContext> toReject = new HashSet<>();

            if (gv.type.equals("NAHR-INS")) {
                for (VariantContext vc : clusteredVariants.get(gv)) {
                    if (!isSameType(gv, vc)) {
                        toReject.add(vc);
                    } else if (gv.seqIndex != getSequenceIndex(vc)) {
                        toReject.add(vc);
                    }
                }

                clusteredVariants.get(gv).removeAll(toReject);
            }
        }

        return filteredVariants;
    }

    private int getSequenceIndex(VariantContext vc) {
        return REFERENCES.get(vc.getAttributeAsString("BACKGROUND", REFERENCES.keySet().iterator().next())).getReferenceSequence().getSequenceDictionary().getSequenceIndex(vc.getContig()) - 1;
    }

    private boolean isSameType(GeneratedVariant gv, VariantContext v) {
        if (gv.type.equalsIgnoreCase("SNV") && v.isSNP() || (v.isIndel() && v.getAlternateAllele(0).length() < 10)) { return true; }

        if (Math.abs(gv.newAllele.length() - gv.oldAllele.length()) <= 100) {
            if (gv.type.equalsIgnoreCase("TD") && v.isSimpleInsertion()) { return true; }
            if (gv.type.equalsIgnoreCase("INS") && v.isSimpleInsertion()) { return true; }
            if (gv.type.equalsIgnoreCase("DEL") && v.isSimpleDeletion()) { return true; }
            if (gv.type.equalsIgnoreCase("STR_EXP") && v.isIndel()) { return true; }
            if (gv.type.equalsIgnoreCase("STR_CON") && v.isIndel()) { return true; }
            if (gv.type.equalsIgnoreCase("INV") && v.isMNP()) { return true; }
            if (gv.type.equalsIgnoreCase("MNP") && v.isMNP()) { return true; }
        } else {
            if (gv.type.equalsIgnoreCase("TD") && (v.isSimpleInsertion() || v.isSymbolic())) { return true; }
            if (gv.type.equalsIgnoreCase("INS") && (v.isSimpleInsertion() || v.isSymbolic())) { return true; }
            if (gv.type.equalsIgnoreCase("DEL") && (v.isSimpleDeletion() || v.isSymbolic())) { return true; }
            if (gv.type.equalsIgnoreCase("STR_EXP") && (v.isIndel() || v.isSymbolic())) { return true; }
            if (gv.type.equalsIgnoreCase("STR_CON") && (v.isIndel() || v.isSymbolic())) { return true; }
            if (gv.type.equalsIgnoreCase("INV") && (v.isMNP() || v.isSymbolic())) { return true; }
            if (gv.type.equalsIgnoreCase("MNP") && (v.isMNP() || v.isSymbolic())) { return true; }
        }

        if (gv.type.equalsIgnoreCase("NAHR-INS") && v.isSymbolic() && v.getAlternateAllele(0).getDisplayString().contains(":")) { return true; }
        if (gv.type.equalsIgnoreCase("NAHR-DEL") && v.isSymbolic()) { return true; }

        return false;
    }

    private Map<GeneratedVariant, Set<VariantContext>> clusterCalledVariants(Set<VariantContext> vcs, List<GeneratedVariant> gvs) {
        Map<GeneratedVariant, Set<VariantContext>> clusteredVariants = new HashMap<>();

        GeneratedVariant reject = new GeneratedVariant("rejected", 0, 0, "", "");
        clusteredVariants.put(reject, new HashSet<>());

        Set<VariantContext> used = new HashSet<>();

        for (GeneratedVariant gv : gvs) {
            clusteredVariants.put(gv, new HashSet<>());

            for (VariantContext v : vcs) {
//                String generatedChildHap = (gv.seedLeft + gv.newAllele + gv.seedRight).replaceAll("\\.", "").toUpperCase();
//                String calledChildHap = (v.getAttributeAsString("CHILD_HAP", "")).toUpperCase();
//
//                if (gv.type.equals("NAHR-INS")) {
//                    generatedChildHap = gv.newAllele.toUpperCase();
//                }
//
//                String calledHap;
//                IndexedReference ref1 = (v.getContig().contains("HB3")) ? REFERENCES.get("HB3") : REFERENCES.get("DD2");
//                if (v.isSymbolic() && v.getAlternateAllele(0).getDisplayString().contains(":")) {
//                    String q = v.getAlternateAllele(0).getDisplayString()
//                            .replaceAll("(HB3|DD2):", "")
//                            .replaceAll("[\\[\\]]", "")
//                            .replaceAll("^[ACGT]", "")
//                            .replaceAll("[ACGT]$", "")
//                            .replaceAll("-.*", "");
//
//                    String[] qa = q.split(":");
//                    String chr2 = qa[0];
//                    long start2 = Integer.parseInt(qa[1]);
//
//                    IndexedReference ref2 = chr2.contains("HB3") ? REFERENCES.get("HB3") : REFERENCES.get("DD2");
//
//                    String calledHap1 = ref1.getReferenceSequence().getSubsequenceAt(v.getContig(), v.getStart() - 100, v.getStart() + 100).getBaseString();
//                    String calledHap2 = ref2.getReferenceSequence().getSubsequenceAt(chr2, start2 - 100, start2 + 100).getBaseString();
//                }
//
//                Set<CanonicalKmer> kmersG = kmerize(generatedChildHap, 21);
//                Set<CanonicalKmer> kmersC = kmerize(calledChildHap, 21);
//                Set<CanonicalKmer> kmersA = new HashSet<>();
//                kmersA.addAll(kmersG);
//                kmersA.addAll(kmersC);
//
//                int numKmersOverlap = 0, numKmers = 0;
//                for (CanonicalKmer kmer : kmersA) {
//                    if (kmersG.contains(kmer)) {
//                        numKmers++;
//
//                        if (kmersC.contains(kmer)) {
//                            numKmersOverlap++;
//                        }
//                    }
//                }
//
//                double pctOverlap = (double) numKmersOverlap / (double) numKmers;

//                String bg = v.getAttributeAsString("BACKGROUND", REFERENCES.keySet().iterator().next());
//                int seqIndex = REFERENCES.get(bg).getReferenceSequence().getSequenceDictionary().getSequenceIndex(v.getContig());

                //if (pctOverlap > 0.50 || isApproxSameLocus(gv, v, seqIndex)) {
                if (isSimilar(gv, v)) {
                    clusteredVariants.get(gv).add(v);
                    used.add(v);
                }
            }
        }

        for (VariantContext v : vcs) {
            if (!used.contains(v)) {
                clusteredVariants.get(reject).add(v);
            }
        }

        return clusteredVariants;
    }

    private boolean isSimilar(GeneratedVariant gv, VariantContext v) {
        boolean sameSeqIndex = false;
        boolean similarPosition = false;
        if (gv.getType().contains("NAHR") && v.getAlternateAllele(0).getDisplayString().contains(":") && !v.getContig().contains("contig") && !v.getAlternateAllele(0).getDisplayString().contains("contig") && !v.getAlternateAllele(0).getDisplayString().contains("TANDEM")) {
            String q = v.getAlternateAllele(0).getDisplayString()
                    .replaceAll("(HB3|DD2):", "")
                    .replaceAll("[\\[\\]]", "")
                    .replaceAll("^[ACGT]+", "")
                    .replaceAll("[ACGT]+$", "")
                    .replaceAll("-.*", "");

            String[] qa = q.split(":");
            String chr2 = qa[0];
            long start2 = Integer.parseInt(qa[1]);

            int nahrSeqIndexGV1 = gv.start.getContig().contains("HB3") ?
                    REFERENCES.get("HB3").getReferenceSequence().getSequence(gv.start.getContig()).getContigIndex() :
                    REFERENCES.get("DD2").getReferenceSequence().getSequence(gv.start.getContig()).getContigIndex();

            int nahrSeqIndexGV2 = gv.stop.getContig().contains("HB3") ?
                    REFERENCES.get("HB3").getReferenceSequence().getSequence(gv.stop.getContig()).getContigIndex() :
                    REFERENCES.get("DD2").getReferenceSequence().getSequence(gv.stop.getContig()).getContigIndex();

            int nahrSeqIndexV1 = v.getContig().contains("HB3") ?
                    REFERENCES.get("HB3").getReferenceSequence().getSequence(v.getContig()).getContigIndex() :
                    REFERENCES.get("DD2").getReferenceSequence().getSequence(v.getContig()).getContigIndex();

            int nahrSeqIndexV2 = chr2.contains("HB3") ?
                    REFERENCES.get("HB3").getReferenceSequence().getSequence(chr2).getContigIndex() :
                    REFERENCES.get("DD2").getReferenceSequence().getSequence(chr2).getContigIndex();

            sameSeqIndex = (nahrSeqIndexGV1 == nahrSeqIndexV1 && nahrSeqIndexGV2 == nahrSeqIndexV2) ||
                              (nahrSeqIndexGV1 == nahrSeqIndexV2 && nahrSeqIndexGV2 == nahrSeqIndexV1);

            similarPosition = (Math.abs(gv.start.getStart() - v.getStart()) < 10000 && Math.abs(gv.start.getStart() - start2) < 10000) ||
                              (Math.abs(gv.stop.getStart() - v.getStart()) < 10000 && Math.abs(gv.stop.getStart() - start2) < 10000);
        } else {
            int seqIndex = -1;
            if (v.getContig().contains("HB3")) {
                seqIndex = REFERENCES.get("HB3").getReferenceSequence().getSequenceDictionary().getSequenceIndex(v.getContig());
            } else {
                seqIndex = REFERENCES.get("DD2").getReferenceSequence().getSequenceDictionary().getSequenceIndex(v.getContig());
            }

            sameSeqIndex = (gv.seqIndex + 1 == seqIndex);

            if (gv.getType().equals("SNV")) {
                similarPosition = Math.abs(gv.start.getEnd() - v.getStart()) < 10 || Math.abs(gv.stop.getStart() - v.getStart()) < 10;
            } else {
                similarPosition = Math.abs(gv.start.getEnd() - v.getStart()) < 500 || Math.abs(gv.stop.getStart() - v.getStart()) < 500;
            }
        }

        boolean similarType = isSameType(gv, v);

        return sameSeqIndex && similarType && similarPosition;
    }

    private Set<CanonicalKmer> kmerize(String seq, int kmerSize) {
        Set<CanonicalKmer> kmers = new HashSet<>();

        for (int i = 0; i <= seq.length() - kmerSize; i++) {
            String sk = seq.substring(i, i + kmerSize).toUpperCase();

            kmers.add(new CanonicalKmer(sk));
        }

        return kmers;
    }

    private List<GeneratedVariant> updateLociOfSimulatedVariants(List<GeneratedVariant> vs) {
        List<GeneratedVariant> nvs = new ArrayList<>();
        for (GeneratedVariant v : vs) {
            String oldHap = (v.seedLeft + v.oldAllele + v.seedRight).toUpperCase();
            if (v.type.equalsIgnoreCase("NAHR-INS")) {
                GeneratedVariant nv1 = new GeneratedVariant(v.type, v.seqIndex, v.start.getStart(), v.oldAllele, v.newAllele, v.seedLeft, v.seedRight, v.start.getContig().contains("HB3") ? "HB3" : "DD2", v.start, v.start, v.posIndex);
                GeneratedVariant nv2 = new GeneratedVariant(v.type, v.seqIndex, v.stop.getStart(), v.oldAllele, v.newAllele, v.seedLeft, v.seedRight, v.stop.getContig().contains("HB3") ? "HB3" : "DD2", v.stop, v.stop, v.posIndex);

                nvs.add(nv1);
                nvs.add(nv2);
            } else {
                IndexedReference ref = REFERENCES.get(v.parent);

                ReferenceSequence rseq = ref.getReferenceSequence().getSequence(ref.getReferenceSequence().getSequenceDictionary().getSequence(v.getSeqIndex() + 1).getSequenceName());
                String seq = rseq.getBaseString().toUpperCase();

                int firstPos = seq.indexOf(oldHap);
                int lastPos = seq.lastIndexOf(oldHap);

                if (firstPos == lastPos && firstPos >= 0) {
                    GeneratedVariant nv = new GeneratedVariant(v.type, v.seqIndex, firstPos + v.seedLeft.length(), v.oldAllele, v.newAllele, v.seedLeft, v.seedRight, v.parent, v.start, v.stop, v.index);
                    v.loci = Arrays.asList(new Pair<>(
                            new Interval(rseq.getName(), firstPos + v.seedLeft.length(), firstPos + v.seedLeft.length()),
                            new Interval(rseq.getName(), firstPos + v.seedLeft.length() + v.oldAllele.length(), firstPos + v.seedLeft.length() + v.oldAllele.length())
                    ));

                    nvs.add(nv);
                }
            }
        }

        return nvs;
    }

    private List<Set<GeneratedVariant>> permuteVariants(Set<GeneratedVariant> vs) {
        List<Set<GeneratedVariant>> gsl = new ArrayList<>();

        List<GeneratedVariant> vl = new ArrayList<>();
        List<Integer> indices = new ArrayList<>();
        int i = 0;
        for (GeneratedVariant v : vs) {
            vl.add(v);
            indices.add(i);
            i++;
        }

        List<List<Integer>> lss = generateCombinatoricLists(indices);
        List<Integer> lsall = new ArrayList<>();
        for (int j = 0; j < indices.size(); j++) {
            lsall.add(j);
        }
        lss.add(0, lsall);

        for (List<Integer> ls : lss) {
            //log.info("bef {}", Joiner.on(", ").join(ls));
        }

        Set<String> alreadyPresent = new HashSet<>();
        Set<Integer> toRemove = new TreeSet<>();
        for (int j = 0; j < lss.size(); j++) {
            List<Integer> ls = lss.get(j);
            String joined = Joiner.on(",").join(ls);

            if (alreadyPresent.contains(joined)) {
                toRemove.add(j);
            }

            alreadyPresent.add(joined);
        }

        List<Integer> toRemoveList = new ArrayList<>(toRemove);
        for (int j = toRemoveList.size() - 1; j >= 0; j--) {
            int ir = toRemoveList.get(j);
            lss.remove(ir);
        }

        for (List<Integer> ls : lss) {
            //log.info("aft {}", Joiner.on(", ").join(ls));
            Set<GeneratedVariant> gs = new HashSet<>();
            for (int l : ls) {
                gs.add(vl.get(l));
            }

            gsl.add(gs);
        }

        return gsl;
    }

    private List<List<Integer>> generateCombinatoricLists(List<Integer> indices) {
        List<List<Integer>> loos = new ArrayList<>();

        for (int i = 0; i < indices.size(); i++) {
            List<Integer> loo = new ArrayList<>();

            for (int j = 0; j < indices.size(); j++) {
                if (i != j) {
                    loo.add(indices.get(j));
                }
            }

            if (loo.size() > 0) {
                loos.add(loo);
            }

            if (loo.size() > 1) {
                loos.addAll(generateCombinatoricLists(loo));
            }
        }

        return loos;
    }

    private List<Pair<List<String>, List<Integer>>> recombine(Map<Integer, List<Map<String, String>>> recombs) {
        /*
        FastaSequenceFile f = new FastaSequenceFile(new File("child.new.fasta"), true);
        ReferenceSequence r = f.nextSequence();

        List<String> rs = new ArrayList<>();
        StringBuilder sbr = new StringBuilder();
        boolean upperCase = isUpperCase(r.getBaseString().substring(0, 1));
        for (int i = 0; i < r.length(); i++) {
            String s = r.getBaseString().substring(i, i + 1);

            if (isUpperCase(s) != upperCase) {
                rs.add(sbr.toString());
                sbr = new StringBuilder();
                upperCase = !upperCase;
            }

            sbr.append(s);
        }

        rs.add(sbr.toString());
        */

        Map<String, Integer> switchIds = new HashMap<>();
        int switchId = 0;
        for (String refName : REFERENCES.keySet()) {
            switchIds.put(refName, switchId);
            switchId++;
        }

        List<Pair<List<String>, List<Integer>>> seqs = new ArrayList<>();

        for (int chrIndex : recombs.keySet()) {
            List<String> seq = new ArrayList<>();
            List<Integer> switches = new ArrayList<>();

            for (Map<String, String> te : recombs.get(chrIndex)) {
                StringBuilder sb = new StringBuilder();

                IndexedFastaSequenceFile ref = REFERENCES.get(te.get("parent")).getReferenceSequence();
                String chr = ref.getSequenceDictionary().getSequence(Integer.valueOf(te.get("chr"))).getSequenceName();
                int start = Integer.valueOf(te.get("start")) + 1;
                int stop = Integer.valueOf(te.get("stop"));

                if (stop > ref.getSequence(chr).length()) {
                    stop = ref.getSequence(chr).length();
                }

                sb.append(ref.getSubsequenceAt(chr, start, stop).getBaseString());

                seq.add(sb.toString());
                switches.add(switchIds.get(te.get("parent")));
            }

            seqs.add(new Pair<>(seq, switches));
        }

        return seqs;
    }

    private List<String> collapseG(List<Pair<List<String>, List<Integer>>> seqs, Set<GeneratedVariant> vs) {
        List<GeneratedVariant> lvs = new ArrayList<>(vs);
        List<StringBuilder> newSbs = new ArrayList<>();
        for (int i = 0; i < seqs.size(); i++) {
            StringBuilder sbs = new StringBuilder();
            for (int j = 0; j < seqs.get(i).getSecond().size(); j++) {
                String sb = seqs.get(i).getFirst().get(j);
                int sw = seqs.get(i).getSecond().get(j);

                sb = sw == 1 ? sb.toUpperCase() : sb.toLowerCase();

                sbs.append(sb);
            }

            newSbs.add(sbs);
        }

        lvs.sort(Comparator.comparingInt(GeneratedVariant::getPosIndex));

        for (int i = lvs.size() - 1; i >= 0; i--) {
            GeneratedVariant gv = lvs.get(i);
            StringBuilder sb = newSbs.get(gv.seqIndex);

            String altHap = gv.seedLeft + gv.oldAllele + gv.seedRight;

            int pos = sb.indexOf(altHap);
            if (pos >= 0) {
                sb = sb.replace(pos + gv.seedLeft.length(), pos + gv.seedLeft.length() + gv.oldAllele.length(), gv.newAllele);
            } else {
                throw new CortexJDKException("WTF?");
            }

            //log.info("");
        }

        List<String> newSeqs = new ArrayList<>();
        for (StringBuilder sb : newSbs) {
            newSeqs.add(sb.toString());
        }

        return newSeqs;
    }

    private List<String> collapseV(Set<VariantContext> vcs) {
        List<VariantContext> lvs = new ArrayList<>(vcs);

        Map<String, StringBuilder> newSbs = new LinkedHashMap<>();
        for (String refName : REFERENCES.keySet()) {
            ReferenceSequence rseq;
            while((rseq = REFERENCES.get(refName).getReferenceSequence().nextSequence()) != null) {
                newSbs.put(rseq.getName(), new StringBuilder(rseq.getBaseString()));
            }
        }

        lvs.sort(Comparator.comparingInt(VariantContext::getStart));

        Map<VariantContext, Integer> newPositions = new HashMap<>();
        for (int i = lvs.size() - 1; i >= 0; i--) {
            VariantContext vc = lvs.get(i);
            StringBuilder sb = newSbs.get(vc.getContig());

            sb = sb.replace(vc.getStart(), vc.getEnd(), vc.getAltAlleleWithHighestAlleleCount().getBaseString());

            for (int j = i + 1; j < lvs.size(); j++) {
                if (lvs.get(j).getContig().equalsIgnoreCase(vc.getContig())) {
                    //log.info("");
                }
            }
        }

        List<String> newSeqs = new ArrayList<>();
        for (String chromName : newSbs.keySet()) {
            StringBuilder sb = newSbs.get(chromName);

            newSeqs.add(sb.toString());
        }

        return newSeqs;
    }
}

