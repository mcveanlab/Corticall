package uk.ac.ox.well.indiana.commands.simulate;

import cern.jet.random.Empirical;
import com.google.common.base.Joiner;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.math3.random.EmpiricalDistribution;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class simchild extends Module {
    @Argument(fullName="reference", shortName="R", doc="Reference sequence")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="vcf", shortName="V", doc="Variants for second parent")
    public VCFFileReader VCF;

    @Argument(fullName="seed", shortName="s", doc="Seed for RNG", required=false)
    public Long SEED;

    @Output
    public File out;

    @Output(fullName="statsOut", shortName="so", doc="The stats output file")
    public PrintStream sout;

    private Random rng;

    Map<String, double[]> empDists;
    Map<String, EmpiricalDistribution> empRates;

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

    @Override
    public void execute() {
        if (SEED == null) { SEED = System.currentTimeMillis(); }

        rng = new Random(SEED);
        initializeRates();

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
            if (!variants.containsKey(vc.getChr())) {
                variants.put(vc.getChr(), new HashMap<Integer, VariantContext>());
            }

            variants.get(vc.getChr()).put(vc.getStart() - 1, vc);

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

        VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder();
        vcwb.setOutputFile(out);
        vcwb.setReferenceDictionary(VCF.getFileHeader().getSequenceDictionary());
        VariantContextWriter vcw = vcwb.build();

        vcw.writeHeader(VCF.getFileHeader());

        for (SAMSequenceRecord ssr : VCF.getFileHeader().getSequenceDictionary().getSequences()) {
            String chr = ssr.getSequenceName();

            if (copiedVariants.containsKey(chr)) {
                for (Integer i : copiedVariants.get(chr).keySet()) {
                    VariantContext vc = copiedVariants.get(chr).get(i);

                    vcw.add(vc);
                }
            }
        }

        vcw.close();
    }
}
