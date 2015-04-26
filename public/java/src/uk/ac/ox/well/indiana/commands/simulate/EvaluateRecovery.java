package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTable;
import uk.ac.ox.well.indiana.utils.containers.DataTables;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.AlignmentUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class EvaluateRecovery extends Module {
    @Argument(fullName="truth", shortName="t", doc="Truth VCF")
    public VCFFileReader TRUTH;

    @Argument(fullName="eval", shortName="e", doc="Eval VCF")
    public VCFFileReader EVAL;

    @Argument(fullName="metrics", shortName="m", doc="Contig metrics")
    public File METRICS;

    @Argument(fullName="bam0", shortName="b0", doc="BAM for parent 0")
    public SAMFileReader BAM0;

    @Argument(fullName="bam1", shortName="b1", doc="BAM for parent 1")
    public SAMFileReader BAM1;

    @Argument(fullName="novelKmerVariantMap", shortName="n", doc="Novel kmer variant map")
    public File NOVEL_KMER_VARIANT_MAP;

    @Output
    public PrintStream out;

    @Output(fullName="weirdout", shortName="wo", doc="Weird out")
    public PrintStream wout;

    private class ContigAlignment {
        public Set<SAMRecord> sr0s = new HashSet<SAMRecord>();
        public Set<SAMRecord> sr1s = new HashSet<SAMRecord>();
    }

    private class ContigStatus {
        public boolean shouldHaveBeenCalled = false;
        public boolean wasCalledAndRetained = false;
        public boolean wasCalledAndRejected = false;

        public String toString() {
            return String.format("%b %b %b", shouldHaveBeenCalled, wasCalledAndRejected, wasCalledAndRetained);
        }
    }

    private class VariantInfo {
        private String variantId;
        private String vclass;
        private String vchr;
        private int vstart;
        private int vstop;

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            VariantInfo that = (VariantInfo) o;

            if (vstart != that.vstart) return false;
            if (vstop != that.vstop) return false;
            if (variantId != null ? !variantId.equals(that.variantId) : that.variantId != null) return false;
            if (vchr != null ? !vchr.equals(that.vchr) : that.vchr != null) return false;
            if (vclass != null ? !vclass.equals(that.vclass) : that.vclass != null) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = variantId != null ? variantId.hashCode() : 0;
            result = 31 * result + (vclass != null ? vclass.hashCode() : 0);
            result = 31 * result + (vchr != null ? vchr.hashCode() : 0);
            result = 31 * result + vstart;
            result = 31 * result + vstop;
            return result;
        }
    }

    private Map<CortexKmer, Set<VariantInfo>> loadNovelKmerVariantMap() {
        TableReader tr = new TableReader(NOVEL_KMER_VARIANT_MAP);

        Map<CortexKmer, Set<VariantInfo>> vis = new HashMap<CortexKmer, Set<VariantInfo>>();
        for (Map<String, String> te : tr) {
            CortexKmer ck = new CortexKmer(te.get("kmer"));

            VariantInfo vi = new VariantInfo();
            vi.variantId = te.get("variantId");
            vi.vclass = te.get("vclass");
            vi.vchr = te.get("vchr");
            vi.vstart = Integer.valueOf(te.get("vstart"));
            vi.vstop = Integer.valueOf(te.get("vstop"));

            if (!vis.containsKey(ck)) {
                vis.put(ck, new HashSet<VariantInfo>());
            }
            vis.get(ck).add(vi);
        }

        return vis;
    }

    private Map<String, ContigAlignment> loadAlignedContigs() {
        Map<String, ContigAlignment> alignedContigs = new HashMap<String, ContigAlignment>();

        for (SAMRecord sr : BAM0) {
            if (!alignedContigs.containsKey(sr.getReadName())) {
                alignedContigs.put(sr.getReadName(), new ContigAlignment());
            }

            alignedContigs.get(sr.getReadName()).sr0s.add(sr);
        }

        for (SAMRecord sr : BAM1) {
            if (!alignedContigs.containsKey(sr.getReadName())) {
                alignedContigs.put(sr.getReadName(), new ContigAlignment());
            }

            alignedContigs.get(sr.getReadName()).sr0s.add(sr);
        }

        return alignedContigs;
    }

    private Map<String, Map<String, String>> loadContigAnnotations() {
        TableReader tr = new TableReader(METRICS);

        Map<String, Map<String, String>> tes = new HashMap<String, Map<String, String>>();
        for (Map<String, String> te : tr) {
            String contigName = te.get("contigName");

            tes.put(contigName, te);
        }

        return tes;
    }

    @Override
    public void execute() {
        DataTables report = new DataTables();

        log.info("Loading debugging info...");
        Map<CortexKmer, Set<VariantInfo>> vis = loadNovelKmerVariantMap();

        log.info("Loading contig alignments...");
        Map<String, ContigAlignment> alignedContigs = loadAlignedContigs();
        Map<String, Map<String, String>> annotatedContigs = loadContigAnnotations();

        log.info("Loading truth dataset...");
        Map<String, IntervalTreeMap<VariantContext>> truth = new HashMap<String, IntervalTreeMap<VariantContext>>();
        Map<String, Set<String>> typeIds = new TreeMap<String, Set<String>>();
        Set<String> ids = new HashSet<String>();
        Set<String> filteredVariants = new HashSet<String>();
        Set<String> contigsWithCompleteVariants = new HashSet<String>();
        Set<String> contigsWithIncompleteVariants = new HashSet<String>();
        Set<String> contigsWithVariants = new HashSet<String>();
        Map<String, Set<VariantContext>> idToVcMap = new HashMap<String, Set<VariantContext>>();
        int recordsLoaded = 0;

        DataTable knownStats = new DataTable("knownStats", "known stats");
        knownStats.addColumns("type", "complete", "partial", "total");
        report.addTable(knownStats);

        for (VariantContext vc : TRUTH) {
            String type = vc.getAttributeAsString("denovo", "unknown");

            String simchildid = vc.getAttributeAsString("id", "unknown");
            if (!idToVcMap.containsKey(simchildid)) {
                idToVcMap.put(simchildid, new HashSet<VariantContext>());
            }
            idToVcMap.get(simchildid).add(vc);

            if (vc.hasAttribute("denovo") && !type.equals("unknown")) {
                Interval interval = new Interval(vc.getChr(), vc.getStart(), vc.getEnd());
                if (!truth.containsKey(vc.getChr())) {
                    truth.put(vc.getChr(), new IntervalTreeMap<VariantContext>());
                }
                truth.get(vc.getChr()).put(interval, vc);

                String id = vc.getAttributeAsString("id", "unknown");
                if (!typeIds.containsKey(type)) {
                    typeIds.put(type, new HashSet<String>());
                }
                typeIds.get(type).add(id);

                if (vc.isNotFiltered()) {
                    ids.add(id);

                    knownStats.set(type, "type", type);
                    knownStats.increment(type, "complete");

                    knownStats.set("all", "type", "all");
                    knownStats.increment("all", "complete");

                    contigsWithCompleteVariants.add(vc.getChr());
                } else {
                    filteredVariants.add(id);
                    knownStats.increment(type, "partial");

                    knownStats.set("all", "type", "all");
                    knownStats.increment("all", "partial");

                    contigsWithIncompleteVariants.add(vc.getChr());
                }

                contigsWithVariants.add(vc.getChr());

                knownStats.increment(type, "total");
                knownStats.increment("all", "total");

                recordsLoaded++;

            }
        }

        DataTable callStats = new DataTable("callStats", "call stats");
        callStats.addColumns("type", "complete", "partial", "rejected", "total");
        report.addTable(callStats);

        for (VariantContext vc : EVAL) {
            if (vc.getAttributeAsBoolean("DENOVO", false)) {
                String event = vc.getAttributeAsString("event", "unknown");

                if (!event.equals("unknown")) {
                    callStats.set(event, "type", event);
                    callStats.increment(event, "total");

                    if (vc.isNotFiltered()) {
                        callStats.increment(event, "complete");
                    } else if (vc.getFilters().contains("PARTIAL")) {
                        callStats.increment(event, "partial");
                    } else {
                        callStats.increment(event, "rejected");
                    }
                }
            }
        }

        Set<String> allIds = new HashSet<String>();
        allIds.addAll(ids);
        allIds.addAll(filteredVariants);

        log.info("  loaded {} records for {} variants ({} more filtered) for {} total (some ids will appear both filtered and unfiltered)", recordsLoaded, ids.size(), filteredVariants.size(), allIds.size());
        log.info("  {} contigs with complete variants, {} with incomplete variants, total of {} contigs with variants", contigsWithCompleteVariants.size(), contigsWithIncompleteVariants.size(), contigsWithVariants.size());

        log.info("Loading eval dataset...");

        Map<String, IntervalTreeMap<VariantContext>> eval = new HashMap<String, IntervalTreeMap<VariantContext>>();

        for (SAMSequenceRecord ssr : EVAL.getFileHeader().getSequenceDictionary().getSequences()) {
            eval.put(ssr.getSequenceName(), new IntervalTreeMap<VariantContext>());
        }

        DataTable recoveryStats = new DataTable("recoveryStats", "recoveryStats");
        report.addTable(recoveryStats);
        recoveryStats.addColumns("type", "tp", "fp", "tn", "fn", "sens", "spec", "fpr", "fdr", "fn1", "fn2", "total1", "total2");
        recoveryStats.set("contig", "type", "contig");
        recoveryStats.set("contig", "tp", 0l);
        recoveryStats.set("contig", "fp", 0l);
        recoveryStats.set("contig", "tn", 0l);
        recoveryStats.set("contig", "fn", 0l);
        recoveryStats.set("contig", "fn1", 0l);
        recoveryStats.set("contig", "fn2", 0l);

        Set<VariantContext> fpSNPs = new HashSet<VariantContext>();
        Set<VariantContext> fpINSs = new HashSet<VariantContext>();

        for (VariantContext vc : EVAL) {
            if (vc.getAttributeAsBoolean("DENOVO", false)) {
                Interval interval = new Interval(vc.getChr(), vc.getStart(), vc.getEnd());
                if (!vc.isSNP()) {
                    interval = new Interval(vc.getChr(), vc.getStart() - 5, vc.getEnd() + 5);
                }
                String event = vc.getAttributeAsString("event", "unknown");

                if (event != null && !event.equals("GC")) {
                    recoveryStats.set(event, "type", event);

                    if (vc.isNotFiltered()) {
                        boolean isTp = false;
                        if (truth.containsKey(vc.getChr()) && truth.get(vc.getChr()).containsOverlapping(interval)) {
                            for (VariantContext knownVariant : truth.get(vc.getChr()).getOverlapping(interval)) {
                                if (knownVariant.isNotFiltered() && knownVariant.getType() == vc.getType()) {
                                    recoveryStats.increment(event, "tp");
                                    isTp = true;

                                    break;
                                }
                            }
                        }

                        if (!isTp) {
                            recoveryStats.increment(event, "fp");

                            if (vc.isSNP()) {
                                fpSNPs.add(vc);
                            } else if (vc.isSimpleInsertion()) {
                                fpINSs.add(vc);
                                wout.println(vc.getChr());
                            }
                        }
                    } else {
                        boolean isFn = false;

                        if (truth.containsKey(vc.getChr()) && truth.get(vc.getChr()).containsOverlapping(interval)) {
                            for (VariantContext knownVariant : truth.get(vc.getChr()).getOverlapping(interval)) {
                                if (knownVariant.isNotFiltered() && knownVariant.getType() == vc.getType()) {
                                    recoveryStats.increment(event, "fn1");
                                    isFn = true;

                                    break;
                                }
                            }
                        }

                        if (!isFn) {
                            recoveryStats.increment(event, "tn");
                        }
                    }

                    recoveryStats.increment(event, "total1");

                    if (!eval.containsKey(vc.getChr())) {
                        eval.put(vc.getChr(), new IntervalTreeMap<VariantContext>());
                    }

                    eval.get(vc.getChr()).put(interval, vc);
                }
            }
        }

        for (IntervalTreeMap<VariantContext> vcMap : truth.values()) {
            for (VariantContext vc : vcMap.values()) {
                String contigName = vc.getChr();
                Interval knownInterval = new Interval(vc.getChr(), vc.getStart(), vc.getEnd());
                String event = vc.getAttributeAsString("denovo", "unknown");

                if (event != null && !event.equals("unknown") && !event.equals("GC")) {
                    if (!eval.containsKey(contigName) || !eval.get(contigName).containsOverlapping(knownInterval)) {
                        recoveryStats.increment(event, "fn2");
                    }

                    recoveryStats.increment(event, "total2");
                }
            }
        }

        Map<String, ContigStatus> contigStatus = new HashMap<String, ContigStatus>();

        for (SAMSequenceRecord ssr : EVAL.getFileHeader().getSequenceDictionary().getSequences()) {
            String contigName = ssr.getSequenceName();

            contigStatus.put(contigName, new ContigStatus());
        }

        for (String contigName : truth.keySet()) {
            boolean atLeastOneIsCallable = false;
            for (VariantContext vc : truth.get(contigName).values()) {
                if (vc.isNotFiltered() && vc.getAttributeAsBoolean("isComplete", false)) {
                    atLeastOneIsCallable = true;
                }
            }

            contigStatus.get(contigName).shouldHaveBeenCalled = atLeastOneIsCallable;
        }

        for (String contigName : eval.keySet()) {
            boolean calledAtLeastOne = eval.get(contigName).size() > 0;
            boolean retainedAtLeastOne = false;

            for (VariantContext vc : eval.get(contigName).values()) {
                if (vc.isNotFiltered() && vc.getAttributeAsBoolean("DENOVO", false)) {
                    retainedAtLeastOne = true;
                }
            }

            if (calledAtLeastOne) {
                if (retainedAtLeastOne) {
                    contigStatus.get(contigName).wasCalledAndRetained = true;
                } else {
                    contigStatus.get(contigName).wasCalledAndRejected = true;
                }
            }
        }

        Set<String> fpContigs = new HashSet<String>();
        Set<String> fn1Contigs = new HashSet<String>();
        Set<String> fn2Contigs = new HashSet<String>();

        for (String contigName : contigStatus.keySet()) {
            ContigStatus cs = contigStatus.get(contigName);

            if (cs.shouldHaveBeenCalled) {
                if (cs.wasCalledAndRetained) {
                    recoveryStats.increment("contig", "tp");
                } else if (cs.wasCalledAndRejected) {
                    recoveryStats.increment("contig", "fn1");
                    fn1Contigs.add(contigName);
                } else {
                    recoveryStats.increment("contig", "fn2");
                    fn2Contigs.add(contigName);
                }

                recoveryStats.increment("contig", "total2");
            } else {
                if (cs.wasCalledAndRetained) {
                    recoveryStats.increment("contig", "fp");
                    fpContigs.add(contigName);
                } else if (cs.wasCalledAndRejected) {
                    recoveryStats.increment("contig", "tn");
                } else {
                    recoveryStats.increment("contig", "tn");
                }
            }

            recoveryStats.increment("contig", "total1");
        }

        for (String pk : recoveryStats.getPrimaryKeys()) {
            float fn = ((Long) recoveryStats.get(pk, "fn1") + (Long) recoveryStats.get(pk, "fn2"));
            recoveryStats.set(pk, "fn", (int) fn);

            float tp = (Long) recoveryStats.get(pk, "tp");
            float fp = (Long) recoveryStats.get(pk, "fp");
            float tn = recoveryStats.has(pk) ? (Long) recoveryStats.get(pk, "tn") : 0l;

            float sens = tp / (tp + fn);
            float spec = tn / (tn + fp);
            float fpr  = fp / (fp + tn);
            float fdr  = fp / (tp + fp);

            recoveryStats.set(pk, "sens", String.format("%.2f", sens));
            recoveryStats.set(pk, "spec", String.format("%.2f", spec));
            recoveryStats.set(pk, "fpr",  String.format("%.2f", fpr));
            recoveryStats.set(pk, "fdr",  String.format("%.2f", fdr));
        }

        DataTable fn2stats = new DataTable("contig_fn2", "contig_fn2");
        fn2stats.addColumn("total");
        report.addTable(fn2stats);

        for (String contigName : fn2Contigs) {
            ContigAlignment ca = alignedContigs.get(contigName);

            boolean isExplained = false;

            if (ca.sr0s.size() > 1) {
                if (ca.sr1s.size() > 1) {
                    fn2stats.increment("stats", "multiple_in_first_multiple_in_second");
                    isExplained = true;
                } else if (ca.sr1s.size() == 0) {
                    fn2stats.increment("stats", "multiple_in_first_absent_in_second");
                    isExplained = true;
                } else {
                    fn2stats.increment("stats", "multiple_in_first_proper_in_second");
                    isExplained = true;
                }
            } else if (ca.sr0s.size() == 0) {
                if (ca.sr1s.size() > 1) {
                    fn2stats.increment("stats", "absent_in_first_multiple_in_second");
                    isExplained = true;
                } else if (ca.sr1s.size() == 0) {
                    fn2stats.increment("stats", "absent_in_first_absent_in_second");
                    isExplained = true;
                } else {
                    fn2stats.increment("stats", "absent_in_first_proper_in_second");
                    isExplained = true;
                }
            } else {
                if (ca.sr1s.size() > 1) {
                    fn2stats.increment("stats", "proper_in_first_multiple_in_second");
                    isExplained = true;
                } else if (ca.sr1s.size() == 0) {
                    SAMRecord sr0 = ca.sr0s.iterator().next();

                    boolean isClipped = false;
                    for (CigarElement ce : sr0.getCigar().getCigarElements()) {
                        if (ce.getOperator().equals(CigarOperator.S) || ce.getOperator().equals(CigarOperator.H)) {
                            isClipped = true;
                        }
                    }

                    if (isClipped) {
                        Map<String, String> te = annotatedContigs.get(contigName);
                        String seq = te.get("seq");

                        List<CigarElement> fwCigar = AlignmentUtils.getForwardCigar(sr0);

                        IntervalTreeMap<Integer> clipRegions = new IntervalTreeMap<Integer>();

                        if (fwCigar.get(0).getOperator().equals(CigarOperator.S) || fwCigar.get(0).getOperator().equals(CigarOperator.H)) {
                            clipRegions.put(new Interval(contigName, 0, fwCigar.get(0).getLength()), 0);
                        }

                        if (fwCigar.get(fwCigar.size() - 1).getOperator().equals(CigarOperator.S) || fwCigar.get(fwCigar.size() - 1).getOperator().equals(CigarOperator.H)) {
                            int length = fwCigar.get(fwCigar.size() - 1).getLength();
                            clipRegions.put(new Interval(contigName, seq.length() - length, seq.length()), 0);
                        }

                        int numVariants = 0;
                        int numClippedVariants = 0;
                        for (VariantContext vc : truth.get(contigName).values()) {
                            Interval interval = new Interval(contigName, vc.getStart(), vc.getStart() + vc.getAlternateAllele(0).length());

                            if (clipRegions.containsOverlapping(interval) || clipRegions.containsContained(interval)) {
                                numClippedVariants++;
                            }

                            numVariants++;
                        }

                        if (numClippedVariants > 0) {
                            if (numClippedVariants == numVariants) {
                                fn2stats.increment("stats", "all_variants_clipped");
                                isExplained = true;
                            } else {
                                fn2stats.increment("stats", "some_variants_clipped");
                                isExplained = true;
                            }
                        }
                    } else {
                        Map<String, String> te = annotatedContigs.get(contigName);
                        String seq = te.get("seq");
                        String ko = te.get("kmerOrigin");

                        for (VariantContext vc : truth.get(contigName).values()) {
                            int novelKmersContained = 0, novelKmersLeft = 0, novelKmersRight = 0;

                            for (int i = vc.getStart(); i < vc.getEnd(); i++) {
                                if (i < ko.length() && ko.charAt(i) == '.') {
                                    novelKmersContained++;
                                }
                            }

                            for (int i = vc.getStart() - 1; i >= 0; i--) {
                                if (ko.charAt(i) == '.') {
                                    novelKmersLeft++;
                                } else {
                                    break;
                                }
                            }

                            for (int i = vc.getEnd(); i < seq.length(); i++) {
                                if (i < ko.length() && ko.charAt(i) == '.') {
                                    novelKmersRight++;
                                } else {
                                    break;
                                }
                            }

                            if (novelKmersLeft == 0 && novelKmersRight == 0) {
                                fn2stats.increment("stats", "missing_novel_kmers");
                                isExplained = true;
                            } else if (novelKmersContained == 0) {
                                fn2stats.increment("stats", "missing_novel_kmers");
                                isExplained = true;
                            }
                        }
                    }
                } else {
                    fn2stats.increment("stats", "proper_in_first_proper_in_second");
                    isExplained = true;
                }
            }

            fn2stats.increment("stats", "total");

            if (!isExplained) {
                fn2stats.increment("stats", "else");
            }
        }

        DataTable fn1stats = new DataTable("contig_fn1", "contig_fn1");
        fn1stats.addColumn("total");
        report.addTable(fn1stats);

        for (String contigName : fn1Contigs) {
            boolean proximity = false, partial = false;

            for (VariantContext vc : eval.get(contigName).values()) {
                if (vc.isFiltered()) {
                    for (String filter : vc.getFilters()) {
                        if (filter.equals("PROXIMITY")) {
                            proximity = true;
                        }

                        if (filter.equals("PARTIAL")) {
                            partial = true;
                        }
                    }
                }
            }

            if (proximity) {
                if (partial) {
                    fn1stats.increment("stats", "proximity_and_partial");
                } else {
                    fn1stats.increment("stats", "proximity");
                }
            } else if (partial) {
                fn1stats.increment("stats", "partial");
            } else {
                fn1stats.increment("stats", "none");
            }

            fn1stats.increment("stats", "total");
        }

        DataTable fpstats = new DataTable("contig_fp", "contig_fp");
        fpstats.addColumn("total");
        report.addTable(fpstats);

        Map<CortexKmer, Integer> kmerMultiplicity = new HashMap<CortexKmer, Integer>();
        TableReader tr2 = new TableReader("test.txt", new String[] { "kmer", "cov" });
        for (Map<String, String> te2 : tr2) {
            CortexKmer ck = new CortexKmer(te2.get("kmer"));
            int cov = Integer.valueOf(te2.get("cov"));
            kmerMultiplicity.put(ck, cov);
        }

        for (String contigName : fpContigs) {
            Map<String, String> te = annotatedContigs.get(contigName);

            int refNovel = Integer.valueOf(te.get("refNovel"));

            if (refNovel > 0) {
                String seq = te.get("seq");
                String ko = te.get("kmerOrigin");
                int kmerSize = seq.length() - ko.length() + 1;

                int numNovelKmers = 0, numNonUniqueKmers = 0;

                Set<String> variantIds = new HashSet<String>();
                for (int i = 0; i < ko.length(); i++) {
                    if (ko.charAt(i) == '.') {
                        CortexKmer ck = new CortexKmer(seq.substring(i, i + kmerSize));
                        if (vis.containsKey(ck)) {
                            for (VariantInfo vi : vis.get(ck)) {
                                variantIds.add(vi.variantId);
                            }
                        }

                        numNovelKmers++;
                        if (kmerMultiplicity.containsKey(ck) && kmerMultiplicity.get(ck) > 1) {
                            numNonUniqueKmers++;
                        }
                    }
                }

                boolean hasPartialVariants = false;
                if (variantIds.size() > 0) {
                    for (String id : variantIds) {
                        Set<VariantContext> vcs = idToVcMap.get(id);

                        for (VariantContext vc : vcs) {
                            if (vc.getFilters().contains("PARTIAL")) {
                                hasPartialVariants = true;
                            }
                        }
                    }
                }

                if (hasPartialVariants) {
                    fpstats.increment("stats", "hasPartialVariants");
                } else {
                    fpstats.increment("stats", "noPartialVariants");
                }
            } else {
                fpstats.increment("stats", "else");
            }

            fpstats.increment("stats", "total");
        }

        DataTable snpFpStats = new DataTable("snp_fp", "snp_fp");
        snpFpStats.addColumn("total");
        report.addTable(snpFpStats);

        for (VariantContext vc : fpSNPs) {
            int clength = annotatedContigs.get(vc.getChr()).get("seq").length();

            if (vc.getStart() < 5 || vc.getEnd() > clength - 5) {
                snpFpStats.increment("stats", "edge_of_contig");
            } else {
                snpFpStats.increment("stats", "else");
            }

            snpFpStats.increment("stats", "total");
        }

        DataTable insFpStats = new DataTable("ins_fp", "ins_fp");
        insFpStats.addColumn("total");
        report.addTable(insFpStats);

        for (int i = 0; i < 65; i++) {
            insFpStats.addColumn("a_" + i);
        }

        for (VariantContext vc : fpINSs) {
            boolean isAlignmentIssue = false;
            int distanceToKnownVariant = Integer.MAX_VALUE;
            if (truth.containsKey(vc.getChr())) {
                for (VariantContext t : truth.get(vc.getChr()).values()) {
                    int distance = Math.abs(vc.getStart() - t.getStart());

                    if (distance < distanceToKnownVariant) {
                        distanceToKnownVariant = distance;
                    }
                }
            }


            insFpStats.increment("stats", "a_" + distanceToKnownVariant);
            insFpStats.increment("stats", "total");
        }

        report.write(out);
    }
}
