package uk.ac.ox.well.indiana.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.apache.commons.math3.random.EmpiricalDistribution;
import org.apache.ivy.util.StringUtils;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class SimChild extends Module {
    @Argument(fullName="reference", shortName="R", doc="Reference sequence")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="vcf", shortName="V", doc="Variants for second parent")
    public VCFFileReader VCF;

    @Argument(fullName="seed", shortName="s", doc="Seed for RNG", required=false)
    public Long SEED;

    @Argument(fullName="strMap", shortName="m", doc="STR map")
    public File STRS;

    @Output
    public File out;

    @Output(fullName="statsOut", shortName="so", doc="The stats output file")
    public PrintStream sout;

    private Random rng;

    private Map<String, double[]> empDists;
    private Map<String, EmpiricalDistribution> empRates;
    private IntervalTreeMap<String> mask = new IntervalTreeMap<String>();

    private Map<Integer, List<Map<String, String>>> loadStrMap() {
        Map<Integer, List<Map<String, String>>> strs = new TreeMap<Integer, List<Map<String, String>>>();
        TableReader tr = new TableReader(STRS);

        for (Map<String, String> te : tr) {
            Integer strLength = Integer.valueOf(te.get("strLength"));

            if (!strs.containsKey(strLength)) {
                strs.put(strLength, new ArrayList<Map<String, String>>());
            }

            strs.get(strLength).add(te);
        }

        return strs;
    }

    private EmpiricalDistribution loadDistribution(double[] rates) {
        EmpiricalDistribution ed = new EmpiricalDistribution();
        ed.load(rates);
        ed.reSeed(SEED);

        return ed;
    }

    private void initializeRates() {
        empDists = new HashMap<String, double[]>();

        empDists.put("Pf3D7_01_v3",	new double[] { 0.842,0.105,0.053 });
        empDists.put("Pf3D7_02_v3",	new double[] { 0.656,0.250,0.031,0.062 });
        empDists.put("Pf3D7_03_v3",	new double[] { 0.415,0.488,0.073,0.024 });
        empDists.put("Pf3D7_04_v3",	new double[] { 0.634,0.366 });
        empDists.put("Pf3D7_05_v3",	new double[] { 0.673,0.269,0.058 });
        empDists.put("Pf3D7_06_v3",	new double[] { 0.511,0.333,0.133,0.022 });
        empDists.put("Pf3D7_07_v3",	new double[] { 0.694,0.163,0.122,0.020 });
        empDists.put("Pf3D7_08_v3",	new double[] { 0.674,0.163,0.140,0.023 });
        empDists.put("Pf3D7_09_v3",	new double[] { 0.600,0.218,0.164,0.018 });
        empDists.put("Pf3D7_10_v3",	new double[] { 0.500,0.385,0.096,0.019 });
        empDists.put("Pf3D7_11_v3",	new double[] { 0.370,0.370,0.167,0.074,0.019 });
        empDists.put("Pf3D7_12_v3",	new double[] { 0.509,0.298,0.070,0.105,0.018 });
        empDists.put("Pf3D7_13_v3",	new double[] { 0.234,0.406,0.297,0.016,0.031,0.016 });
        empDists.put("Pf3D7_14_v3",	new double[] { 0.313,0.254,0.239,0.060,0.090,0.030,0.015 });

        empRates = new HashMap<String, EmpiricalDistribution>();

        for (String chr : empDists.keySet()) {
            empRates.put(chr, loadDistribution(empDists.get(chr)));
        }
    }

    private boolean canRecombine(String chr) { return empDists.containsKey(chr); }

    private boolean shouldRecombine() { return rng.nextBoolean(); }

    private int numRecombs(String chr) {
        if (empRates.containsKey(chr)) {
            double val = empRates.get(chr).getNextValue();

            for (int i = 0; i < empDists.get(chr).length; i++) {
                if (Math.abs(val - empDists.get(chr)[i]) < 0.001) {
                    return i+1;
                }
            }
        }

        return 0;
    }

    private int whereRecombine(int length) { return rng.nextInt(length); }

    private SAMSequenceRecord getRandomAutosome() {
        SAMSequenceRecord ssr;

        do {
            int seqIndex = rng.nextInt(REF.getSequenceDictionary().size());
            ssr = REF.getSequenceDictionary().getSequence(seqIndex);
        } while (ssr.getSequenceName().equals("M76611") || ssr.getSequenceName().equals("PFC10_API_IRAB"));

        return ssr;
    }

    private Allele getRefAllele(String chr, int pos, int length) {
        return Allele.create(new String(REF.getSubsequenceAt(chr, pos, pos + length).getBases()), true);
    }

    private Allele getRandomAltAllele(Allele refAllele, int length) {
        Allele altAllele = Allele.create(refAllele.getBases(), false);

        final String[] bases = {"A", "C", "G", "T"};

        if (length == 0) {
            while (altAllele.basesMatch(refAllele)) {
                int baseIndex = rng.nextInt(bases.length);
                altAllele = Allele.create(bases[baseIndex], false);
            }
        } else {
            StringBuilder altAlleleString = new StringBuilder();
            altAlleleString.append(refAllele.getBaseString());

            for (int i = 0; i < length; i++) {
                int baseIndex = rng.nextInt(bases.length);

                altAlleleString.append(bases[baseIndex]);
            }

            altAllele = Allele.create(altAlleleString.toString(), false);
        }

        return altAllele;
    }

    private void addToMask(String chr, int pos, int length) {
        Interval interval = new Interval(chr, pos - 2, pos + length + 2);

        mask.put(interval, null);
    }

    private boolean isInMask(String chr, int pos, int length) {
        Interval shortInterval = new Interval(chr, pos, pos);
        Interval fullInterval  = new Interval(chr, pos, pos + length);

        return mask.containsOverlapping(shortInterval) || mask.containsOverlapping(fullInterval);
    }

    private void addDeNovoSNPs(Map<String, Map<Integer, List<VariantContext>>> vcs, int numVariants, String sampleName) {
        for (int i = 0; i < numVariants; i++) {
            SAMSequenceRecord ssr = getRandomAutosome();
            int pos;
            do {
                pos = rng.nextInt(ssr.getSequenceLength()) + 1;
            } while (isInMask(ssr.getSequenceName(), pos, 1));

            Allele refAllele = getRefAllele(ssr.getSequenceName(), pos, 0);
            Allele altAllele = getRandomAltAllele(refAllele, 0);

            Genotype newg = (new GenotypeBuilder(sampleName, Arrays.asList(altAllele))).make();

            VariantContext vcn = (new VariantContextBuilder())
                    .chr(ssr.getSequenceName())
                    .start(pos)
                    .stop(pos)
                    .noID()
                    .attribute("DENOVO", "SNP")
                    .alleles(Arrays.asList(refAllele, altAllele))
                    .genotypes(newg)
                    .make();

            if (!vcs.containsKey(ssr.getSequenceName())) {
                vcs.put(ssr.getSequenceName(), new TreeMap<Integer, List<VariantContext>>());
            }

            if (!vcs.get(ssr.getSequenceName()).containsKey(pos)) {
                vcs.get(ssr.getSequenceName()).put(pos, new ArrayList<VariantContext>());
            }

            vcs.get(ssr.getSequenceName()).get(pos).add(vcn);

            addToMask(ssr.getSequenceName(), pos, 1);
        }
    }

    private void addDeNovoInsertions(Map<String, Map<Integer, List<VariantContext>>> vcs, int numVariants, String sampleName, int length) {
        for (int i = 0; i < numVariants; i++) {
            SAMSequenceRecord ssr = getRandomAutosome();
            int pos;
            do {
                pos = rng.nextInt(ssr.getSequenceLength() - length) + 1;
            } while (isInMask(ssr.getSequenceName(), pos, 1));

            Allele refAllele = getRefAllele(ssr.getSequenceName(), pos, 0);
            Allele altAllele;

            String refString = new String(REF.getSubsequenceAt(ssr.getSequenceName(), pos, pos + length).getBases());
            do {
                altAllele = getRandomAltAllele(refAllele, length);
            } while (altAllele.basesMatch(refString));

            Genotype newg = (new GenotypeBuilder(sampleName, Arrays.asList(altAllele))).make();

            VariantContext vcn = (new VariantContextBuilder())
                    .chr(ssr.getSequenceName())
                    .start(pos)
                    .stop(pos)
                    .noID()
                    .attribute("DENOVO", "INS")
                    .alleles(Arrays.asList(refAllele, altAllele))
                    .genotypes(newg)
                    .make();

            if (!vcs.containsKey(ssr.getSequenceName())) {
                vcs.put(ssr.getSequenceName(), new TreeMap<Integer, List<VariantContext>>());
            }

            if (!vcs.get(ssr.getSequenceName()).containsKey(pos)) {
                vcs.get(ssr.getSequenceName()).put(pos, new ArrayList<VariantContext>());
            }

            vcs.get(ssr.getSequenceName()).get(pos).add(vcn);

            addToMask(ssr.getSequenceName(), pos, 1);
        }
    }

    private void addDeNovoDeletions(Map<String, Map<Integer, List<VariantContext>>> vcs, int numVariants, String sampleName, int length) {
        for (int i = 0; i < numVariants; i++) {
            SAMSequenceRecord ssr = getRandomAutosome();
            int pos;
            do {
                pos = rng.nextInt(ssr.getSequenceLength() - length) + 1;
            } while (isInMask(ssr.getSequenceName(), pos, length));

            Allele refAllele = Allele.create(getRefAllele(ssr.getSequenceName(), pos, length).getBases(), true);
            Allele altAllele = Allele.create(refAllele.getBases()[0], false);

            Genotype newg = (new GenotypeBuilder(sampleName, Arrays.asList(altAllele))).make();

            VariantContext vcn = (new VariantContextBuilder())
                    .chr(ssr.getSequenceName())
                    .start(pos)
                    .computeEndFromAlleles(Arrays.asList(refAllele, altAllele), pos)
                    .noID()
                    .attribute("DENOVO", "DEL")
                    .alleles(Arrays.asList(refAllele, altAllele))
                    .genotypes(newg)
                    .make();

            if (!vcs.containsKey(ssr.getSequenceName())) {
                vcs.put(ssr.getSequenceName(), new TreeMap<Integer, List<VariantContext>>());
            }

            if (!vcs.get(ssr.getSequenceName()).containsKey(pos)) {
                vcs.get(ssr.getSequenceName()).put(pos, new ArrayList<VariantContext>());
            }

            vcs.get(ssr.getSequenceName()).get(pos).add(vcn);

            addToMask(ssr.getSequenceName(), pos, length);
        }
    }

    private void addDeNovoInversions(Map<String, Map<Integer, List<VariantContext>>> vcs, int numVariants, String sampleName, int length) {
        for (int i = 0; i < numVariants; i++) {
            SAMSequenceRecord ssr;
            int pos;

            Allele refAllele, altAllele;

            do {
                ssr = getRandomAutosome();
                pos = rng.nextInt(ssr.getSequenceLength() - length) + 1;

                refAllele = getRefAllele(ssr.getSequenceName(), pos, length);
                altAllele = Allele.create(SequenceUtils.reverse(refAllele.getBases()), false);
            } while (refAllele.basesMatch(altAllele) || isInMask(ssr.getSequenceName(), pos, length));

            Genotype newg = (new GenotypeBuilder(sampleName, Arrays.asList(altAllele))).make();

            VariantContext vcn = (new VariantContextBuilder())
                    .chr(ssr.getSequenceName())
                    .start(pos)
                    .computeEndFromAlleles(Arrays.asList(refAllele, altAllele), pos)
                    .noID()
                    .attribute("DENOVO", "INV")
                    .alleles(Arrays.asList(refAllele, altAllele))
                    .genotypes(newg)
                    .make();

            if (!vcs.containsKey(ssr.getSequenceName())) {
                vcs.put(ssr.getSequenceName(), new TreeMap<Integer, List<VariantContext>>());
            }

            if (!vcs.get(ssr.getSequenceName()).containsKey(pos)) {
                vcs.get(ssr.getSequenceName()).put(pos, new ArrayList<VariantContext>());
            }

            vcs.get(ssr.getSequenceName()).get(pos).add(vcn);

            addToMask(ssr.getSequenceName(), pos, length);
        }
    }

    private void addDeNovoStrExpansions(Map<String, Map<Integer, List<VariantContext>>> vcs, int numVariants, String sampleName, int length, Map<Integer, List<Map<String, String>>> strMap) {
        List<Map<String, String>> strs = strMap.get(length);

        for (int i = 0; i < numVariants; i++) {
            String chr;
            int pos;

            Allele refAllele, altAllele;
            String repUnit = null;
            int repsBefore = 0;
            int repsAfter = 0;

            do {
                int strIndex = rng.nextInt(strs.size()) + 1;
                Map<String, String> te = strs.get(strIndex);

                chr = te.get("chr");
                int start = Integer.valueOf(te.get("start"));
                pos = start - 1;
                int numRepeats = rng.nextInt(Integer.valueOf(te.get("numRepeats")) - 1) + 1;

                refAllele = getRefAllele(chr, start - 1, 0);
                altAllele = Allele.create(StringUtils.repeat(te.get("str"), numRepeats), false);

                repUnit = te.get("str");
                repsBefore = Integer.valueOf(te.get("numRepeats"));
                repsAfter = repsBefore + numRepeats;
            } while (isInMask(chr, pos, altAllele.length()));

            //log.info("  {}:{} {} {}", chr, pos, refAllele, altAllele);

            Genotype newg = (new GenotypeBuilder(sampleName, Arrays.asList(altAllele))).make();

            VariantContext vcn = (new VariantContextBuilder())
                    .chr(chr)
                    .start(pos)
                    .computeEndFromAlleles(Arrays.asList(refAllele, altAllele), pos)
                    .noID()
                    .attribute("DENOVO", "STR_EXP")
                    .attribute("REPUNIT", repUnit)
                    .attribute("REPS_BEFORE", repsBefore)
                    .attribute("REPS_AFTER", repsAfter)
                    .alleles(Arrays.asList(refAllele, altAllele))
                    .genotypes(newg)
                    .make();

            if (!vcs.containsKey(chr)) {
                vcs.put(chr, new TreeMap<Integer, List<VariantContext>>());
            }

            if (!vcs.get(chr).containsKey(pos)) {
                vcs.get(chr).put(pos, new ArrayList<VariantContext>());
            }

            vcs.get(chr).get(pos).add(vcn);

            addToMask(chr, pos, altAllele.length());
        }
    }

    private void addDeNovoStrContractions(Map<String, Map<Integer, List<VariantContext>>> vcs, int numVariants, String sampleName, int length, Map<Integer, List<Map<String, String>>> strMap) {
        List<Map<String, String>> strs = strMap.get(length);

        for (int i = 0; i < numVariants; i++) {
            String chr;
            int pos;

            Allele refAllele, altAllele;
            String repUnit = null;
            int repsBefore = 0;
            int repsAfter = 0;

            do {
                int strIndex = rng.nextInt(strs.size()) + 1;
                Map<String, String> te = strs.get(strIndex);

                chr = te.get("chr");
                int start = Integer.valueOf(te.get("start"));
                pos = start - 1;
                int numRepeats = rng.nextInt(Integer.valueOf(te.get("numRepeats")) - 1) + 1;

                String prevBase = getRefAllele(chr, start - 1, 0).getBaseString();

                altAllele = Allele.create(prevBase, false);
                refAllele = Allele.create(prevBase + StringUtils.repeat(te.get("str"), numRepeats), true);

                repUnit = te.get("str");
                repsBefore = Integer.valueOf(te.get("numRepeats"));
                repsAfter = repsBefore - numRepeats;
            } while (isInMask(chr, pos, altAllele.length()));

            //log.info("  {}:{} {} {}", chr, pos, refAllele, altAllele);

            Genotype newg = (new GenotypeBuilder(sampleName, Arrays.asList(altAllele))).make();

            VariantContext vcn = (new VariantContextBuilder())
                    .chr(chr)
                    .start(pos)
                    .computeEndFromAlleles(Arrays.asList(refAllele, altAllele), pos)
                    .noID()
                    .attribute("DENOVO", "STR_CON")
                    .attribute("REPUNIT", repUnit)
                    .attribute("REPS_BEFORE", repsBefore)
                    .attribute("REPS_AFTER", repsAfter)
                    .alleles(Arrays.asList(refAllele, altAllele))
                    .genotypes(newg)
                    .make();

            if (!vcs.containsKey(chr)) {
                vcs.put(chr, new TreeMap<Integer, List<VariantContext>>());
            }

            if (!vcs.get(chr).containsKey(pos)) {
                vcs.get(chr).put(pos, new ArrayList<VariantContext>());
            }

            vcs.get(chr).get(pos).add(vcn);

            addToMask(chr, pos, altAllele.length());
        }
    }

    private void addDeNovoTandemDuplications(Map<String, Map<Integer, List<VariantContext>>> vcs, int numVariants, String sampleName, int length) {
        for (int i = 0; i < numVariants; i++) {
            SAMSequenceRecord ssr;
            int pos;

            Allele refAllele, altAllele;

            do {
                ssr = getRandomAutosome();
                pos = rng.nextInt(ssr.getSequenceLength() - length) + 1;

                refAllele = getRefAllele(ssr.getSequenceName(), pos, 0);
                altAllele = Allele.create(getRefAllele(ssr.getSequenceName(), pos, length + 1).getBases(), false);
            } while (isInMask(ssr.getSequenceName(), pos, length));

            Genotype newg = (new GenotypeBuilder(sampleName, Arrays.asList(altAllele))).make();

            VariantContext vcn = (new VariantContextBuilder())
                    .chr(ssr.getSequenceName())
                    .start(pos)
                    .computeEndFromAlleles(Arrays.asList(refAllele, altAllele), pos)
                    .noID()
                    .attribute("DENOVO", "TD")
                    .alleles(Arrays.asList(refAllele, altAllele))
                    .genotypes(newg)
                    .make();

            if (!vcs.containsKey(ssr.getSequenceName())) {
                vcs.put(ssr.getSequenceName(), new TreeMap<Integer, List<VariantContext>>());
            }

            if (!vcs.get(ssr.getSequenceName()).containsKey(pos)) {
                vcs.get(ssr.getSequenceName()).put(pos, new ArrayList<VariantContext>());
            }

            vcs.get(ssr.getSequenceName()).get(pos).add(vcn);

            addToMask(ssr.getSequenceName(), pos, length);
        }
    }


    @Override
    public void execute() {
        if (SEED == null) { SEED = System.currentTimeMillis(); }

        rng = new Random(SEED);
        initializeRates();

        log.info("Loading STR positions...");
        Map<Integer, List<Map<String, String>>> strMap = loadStrMap();
        for (Integer strLength : strMap.keySet()) {
            log.info("  strLength={}: instances={}", strLength, strMap.get(strLength).size());
        }

        Map<String, SAMSequenceRecord> chrs = new TreeMap<String, SAMSequenceRecord>();
        for (SAMSequenceRecord ssr : REF.getSequenceDictionary().getSequences()) {
            chrs.put(ssr.getSequenceName(), ssr);
        }

        for (String chr : chrs.keySet()) {
            SAMSequenceRecord ssr = chrs.get(chr);

            List<Integer> recombNums = new ArrayList<Integer>();
            for (int i = 0; i < 10000; i++) {
                recombNums.add(numRecombs(ssr.getSequenceName()));
            }

            sout.println(chr + " " + Joiner.on(",").join(recombNums));
        }

        log.info("Simulating recombination breakpoints (seed={})...", SEED);
        Map<String, Set<Integer>> recombs = new TreeMap<String, Set<Integer>>();
        for (String chr : chrs.keySet()) {
            SAMSequenceRecord ssr = chrs.get(chr);

            if (canRecombine(chr) && shouldRecombine()) {
                int recombNum = numRecombs(ssr.getSequenceName());

                Set<Integer> recombPositions = new TreeSet<Integer>();
                for (int i = 0; i < recombNum; i++) {
                    int recombPos = whereRecombine(ssr.getSequenceLength());

                    recombPositions.add(recombPos);
                }

                recombs.put(chr, recombPositions);

                log.info("  {}: length={} numRecombs={} recombs={}", ssr.getSequenceName(), ssr.getSequenceLength(), recombNum, recombPositions);
            } else {
                log.info("  {}: length={} numRecombs=0", ssr.getSequenceName(), ssr.getSequenceLength());
            }
        }

        log.info("Loading variants...");
        int numVariants = 0;
        Map<String, Map<Integer, VariantContext>> variants = new HashMap<String, Map<Integer, VariantContext>>();
        for (VariantContext vc : VCF) {
            VariantContext newvc = (new VariantContextBuilder(vc))
                    .attributes(null)
                    .make();

            if (!variants.containsKey(newvc.getChr())) {
                variants.put(newvc.getChr(), new HashMap<Integer, VariantContext>());
            }

            variants.get(newvc.getChr()).put(newvc.getStart() - 1, newvc);

            numVariants++;
        }
        log.info("  loaded {} variants", numVariants);

        log.info("Constructing base genome...");
        int numCopiedVariants = 0;
        Map<String, Map<Integer, VariantContext>> copiedVariants = new TreeMap<String, Map<Integer, VariantContext>>();
        for (String chr : variants.keySet()) {
            boolean copyFromAlt = rng.nextBoolean();
            for (int i = 0; i < chrs.get(chr).getSequenceLength(); i++) {
                if (recombs.containsKey(chr) && recombs.get(chr).contains(i)) {
                    copyFromAlt = !copyFromAlt;

                    VariantContext recombvc = (new VariantContextBuilder())
                            .chr(chr)
                            .start(i)
                            .stop(i)
                            .alleles(Arrays.asList(getRefAllele(chr, i, 0)))
                            .genotypes(
                                    (new GenotypeBuilder()
                                            .name("child_" + SEED)
                                            .alleles(Arrays.asList(getRefAllele(chr, i, 0))))
                                    .make()
                            )
                            .attribute("DENOVO", "RECOMB")
                            .make();

                    if (!copiedVariants.containsKey(chr)) {
                        copiedVariants.put(chr, new TreeMap<Integer, VariantContext>());
                    }

                    copiedVariants.get(chr).put(i, recombvc);
                }

                if (copyFromAlt && variants.containsKey(chr) && variants.get(chr).containsKey(i)) {
                    if (!copiedVariants.containsKey(chr)) {
                        copiedVariants.put(chr, new TreeMap<Integer, VariantContext>());
                    }

                    copiedVariants.get(chr).put(i, variants.get(chr).get(i));

                    numCopiedVariants++;
                }
            }
        }
        log.info("  copied {} variants", numCopiedVariants);

        log.info("Simulating gene conversion events...");
        int removedVariants = 0, addedVariants = 0;
        for (String chr : variants.keySet()) {
            int numGCEvents = 3;

            List<Integer> positions = new ArrayList<Integer>();
            for (Integer i : variants.get(chr).keySet()) {
                positions.add(i);
            }

            Collections.sort(positions);

            for (int i = 0; i < numGCEvents; i++) {
                int numVariantsInGC = rng.nextInt(2) + 1;

                int index;
                do {
                    index = rng.nextInt(positions.size() - numVariantsInGC - 1);
                } while (index <= 0);

                List<Integer> gcVariants = new ArrayList<Integer>();
                gcVariants.add(positions.get(index));

                for (int j = 1; j < numVariantsInGC; j++) {
                    gcVariants.add(positions.get(index + j));
                }

                for (int gcVariant : gcVariants) {
                    if (copiedVariants.containsKey(chr) && copiedVariants.get(chr).containsKey(gcVariant)) {
                        VariantContext vc = copiedVariants.get(chr).get(gcVariant);
                        VariantContext newvc = (new VariantContextBuilder(vc))
                                .attribute("DENOVO", "GC")
                                .attribute("GCINDEX", i)
                                .alleles(Arrays.asList(vc.getReference()))
                                .genotypes((new GenotypeBuilder(vc.getGenotype(0)))
                                        .alleles(Arrays.asList(vc.getReference()))
                                        .make())
                                .make();

                        //copiedVariants.get(chr).remove(gcVariant);
                        copiedVariants.get(chr).put(gcVariant, newvc);

                        removedVariants++;
                    } else {
                        if (!copiedVariants.containsKey(chr)) {
                            copiedVariants.put(chr, new TreeMap<Integer, VariantContext>());
                        }

                        VariantContext vc = variants.get(chr).get(gcVariant);
                        VariantContext newvc = (new VariantContextBuilder(vc))
                                .attribute("DENOVO", "GC")
                                .attribute("GCINDEX", i)
                                .make();

                        //copiedVariants.get(chr).put(gcVariant, variants.get(chr).get(gcVariant));
                        copiedVariants.get(chr).put(gcVariant, newvc);

                        addedVariants++;
                    }
                }
            }
        }
        log.info("  removed {} variants, added {} variants", removedVariants, addedVariants);

        VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder();
        vcwb.setOutputFile(out);
        vcwb.setReferenceDictionary(VCF.getFileHeader().getSequenceDictionary());
        VariantContextWriter vcw = vcwb.build();

        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.add(VCF.getFileHeader().getFormatHeaderLine("GT"));

        Set<String> sampleNames = new HashSet<String>();
        sampleNames.add("child_" + SEED);

        VCFHeader header = new VCFHeader(headerLines, sampleNames);
        header.setSequenceDictionary(VCF.getFileHeader().getSequenceDictionary());
        header.addMetaDataLine(new VCFInfoHeaderLine("DENOVO", 1, VCFHeaderLineType.String, "The type of de novo event added"));
        header.addMetaDataLine(new VCFInfoHeaderLine("GCINDEX", 1, VCFHeaderLineType.Integer, "The id of the gene conversion event"));
        header.addMetaDataLine(new VCFInfoHeaderLine("REPUNIT", 1, VCFHeaderLineType.String, "The repeated unit in the STR"));
        header.addMetaDataLine(new VCFInfoHeaderLine("REPS_BEFORE", 1, VCFHeaderLineType.Integer, "The number of repeated units in the STR before modification"));
        header.addMetaDataLine(new VCFInfoHeaderLine("REPS_AFTER", 1, VCFHeaderLineType.Integer, "The number of repeated units in the STR after modification"));
        header.addMetaDataLine(new VCFInfoHeaderLine("SIMID", 1, VCFHeaderLineType.String, "The ID of the simulated variant"));
        for (VCFInfoHeaderLine infoHeaderLine : VCF.getFileHeader().getInfoHeaderLines()) {
            header.addMetaDataLine(infoHeaderLine);
        }

        vcw.writeHeader(header);

        Map<String, Map<Integer, List<VariantContext>>> vcs = new HashMap<String, Map<Integer, List<VariantContext>>>();

        for (SAMSequenceRecord ssr : VCF.getFileHeader().getSequenceDictionary().getSequences()) {
            String chr = ssr.getSequenceName();

            if (copiedVariants.containsKey(chr)) {
                for (Integer i : copiedVariants.get(chr).keySet()) {
                    VariantContext vc = copiedVariants.get(chr).get(i);

                    Genotype oldg = vc.getGenotype(0);
                    Genotype newg = (new GenotypeBuilder("child_" + SEED, oldg.getAlleles())).make();
                    GenotypesContext gc = GenotypesContext.create(newg);

                    VariantContext vcn = (new VariantContextBuilder(vc))
                            .noID()
                            .genotypes(gc)
                            .make();

                    if (!vcs.containsKey(ssr.getSequenceName())) {
                        vcs.put(ssr.getSequenceName(), new TreeMap<Integer, List<VariantContext>>());
                    }

                    if (!vcs.get(ssr.getSequenceName()).containsKey(vc.getStart())) {
                        vcs.get(ssr.getSequenceName()).put(vc.getStart(), new ArrayList<VariantContext>());
                    }

                    vcs.get(ssr.getSequenceName()).get(vc.getStart()).add(vcn);
                }
            }
        }

        addDeNovoSNPs(vcs, 250, "child_" + SEED);
        for (int i = 1; i <= 250; i++) {
            addDeNovoInsertions(vcs, 10, "child_" + SEED, i);
            addDeNovoDeletions(vcs,  10, "child_" + SEED, i);
            addDeNovoInversions(vcs, 10, "child_" + SEED, i);
        }

        for (int strLength = 2; strLength <= 5; strLength++) {
            addDeNovoStrExpansions(vcs, 10, "child_" + SEED, strLength, strMap);
            addDeNovoStrContractions(vcs, 10, "child_" + SEED, strLength, strMap);
        }

        for (int tdLength = 10; tdLength <= 50; tdLength++) {
            addDeNovoTandemDuplications(vcs, 10, "child_" + SEED, tdLength);
        }

        int simid = 0;
        for (SAMSequenceRecord ssr : VCF.getFileHeader().getSequenceDictionary().getSequences()) {
            String chr = ssr.getSequenceName();

            if (vcs.containsKey(chr)) {
                for (int pos : vcs.get(chr).keySet()) {
                    for (VariantContext vc : vcs.get(chr).get(pos)) {
                        VariantContext newvc = (new VariantContextBuilder(vc))
                                .attribute("SIMID", "simchild" + simid)
                                .make();

                        vcw.add(newvc);

                        simid++;
                    }
                }
            }
        }

        vcw.close();
    }
}
