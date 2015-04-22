package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.containers.DataTable;

import java.util.*;

public class EvaluateRecovery extends Module {
    @Argument(fullName="truth", shortName="t", doc="Truth VCF")
    public VCFFileReader TRUTH;

    @Argument(fullName="eval", shortName="e", doc="Eval VCF")
    public VCFFileReader EVAL;

    @Override
    public void execute() {
        log.info("Loading truth dataset...");

        Map<String, IntervalTreeMap<VariantContext>> truth = new HashMap<String, IntervalTreeMap<VariantContext>>();
        Map<String, Set<String>> typeIds = new TreeMap<String, Set<String>>();
        Set<String> ids = new HashSet<String>();
        Set<String> filteredVariants = new HashSet<String>();
        Set<String> contigsWithCompleteVariants = new HashSet<String>();
        Set<String> contigsWithIncompleteVariants = new HashSet<String>();
        Set<String> contigsWithVariants = new HashSet<String>();
        int recordsLoaded = 0;

        DataTable knownStats = new DataTable("knownStats", "known stats");
        knownStats.addColumns("type", "complete", "partial", "total");

        for (VariantContext vc : TRUTH) {
            String type = vc.getAttributeAsString("denovo", "unknown");

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

        log.info("\n{}", callStats);


        Set<String> allIds = new HashSet<String>();
        allIds.addAll(ids);
        allIds.addAll(filteredVariants);

        log.info("  loaded {} records for {} variants ({} more filtered) for {} total (some ids will appear both filtered and unfiltered)", recordsLoaded, ids.size(), filteredVariants.size(), allIds.size());
        log.info("  {} contigs with complete variants, {} with incomplete variants, total of {} contigs with variants", contigsWithCompleteVariants.size(), contigsWithIncompleteVariants.size(), contigsWithVariants.size());
        log.info("\n{}", knownStats);

        log.info("Loading eval dataset...");

        Map<String, IntervalTreeMap<VariantContext>> eval = new HashMap<String, IntervalTreeMap<VariantContext>>();

        DataTable recoveryStats = new DataTable("recoveryStats", "recoveryStats");
        recoveryStats.addColumns("type", "tp", "fp", "tn", "fn", "fn1", "fn2", "total1", "total2", "sens", "spec", "fpr", "fdr");
        recoveryStats.set("contig", "type", "contig");
        //recoveryStats.set("GC", "type", "GC");

        Set<String> seenContigs = new HashSet<String>();
        Set<String> fn1Contigs = new HashSet<String>();

        for (VariantContext vc : EVAL) {
            if (vc.getAttributeAsBoolean("DENOVO", false)) {
                if (!seenContigs.contains(vc.getChr())) {
                    if (vc.isNotFiltered() || (vc.getFilters().contains("PARTIAL") && !vc.getFilters().contains("PROXIMITY"))) {
                        if (contigsWithVariants.contains(vc.getChr())) {
                            recoveryStats.increment("contig", "tp");
                        } else {
                            recoveryStats.increment("contig", "fp");
                        }
                    } else {
                        if (!contigsWithVariants.contains(vc.getChr())) {
                            recoveryStats.increment("contig", "tn");
                        } else {
                            recoveryStats.increment("contig", "fn1");

                            fn1Contigs.add(vc.getChr());
                        }
                    }

                    recoveryStats.increment("contig", "total1");

                    seenContigs.add(vc.getChr());
                }

                Interval interval = new Interval(vc.getChr(), vc.getStart(), vc.getEnd());
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

        seenContigs.clear();

        for (String contigName : truth.keySet()) {
            if (!seenContigs.contains(contigName)) {
                recoveryStats.increment("contig", "total2");

                if (!eval.containsKey(contigName) && !fn1Contigs.contains(contigName)) {
                    recoveryStats.increment("contig", "fn2");
                }
            }

            seenContigs.add(contigName);

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

        for (String pk : recoveryStats.getPrimaryKeys()) {
            float fn = ((Long) recoveryStats.get(pk, "fn1") + (Long) recoveryStats.get(pk, "fn2"));
            recoveryStats.set(pk, "fn", (int) fn);

            float tp = (Long) recoveryStats.get(pk, "tp");
            float fp = (Long) recoveryStats.get(pk, "fp");
            float tn = (Long) recoveryStats.get(pk, "tn");

            float sens = tp / (tp + fn);
            float spec = tn / (tn + fp);
            float fpr  = fp / (fp + tn);
            float fdr  = fp / (tp + fp);

            recoveryStats.set(pk, "sens", sens);
            recoveryStats.set(pk, "spec", spec);
            recoveryStats.set(pk, "fpr",  fpr);
            recoveryStats.set(pk, "fdr",  fdr);
        }

        log.info("\n{}", recoveryStats);
        log.info("--");
    }
}
