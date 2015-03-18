package uk.ac.ox.well.indiana.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.ivy.util.StringUtils;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTable;
import uk.ac.ox.well.indiana.utils.containers.DataTables;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.math.MoreMathUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class SimReads extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public FastaSequenceFile REFERENCE;

    @Argument(fullName="readStartDistribution", shortName="c", doc="Read start distribution")
    public File READ_START_DIST;

    @Argument(fullName="errorSimProfile", shortName="e", doc="Error simulation profile")
    public File SIM_PROFILE;

    @Argument(fullName="seed", shortName="s", doc="Seed for RNG", required=false)
    public Long SEED;

    @Output(fullName="fastqEnd1", shortName="o1", doc="Fastq end 1 out")
    public PrintStream o1;

    @Output(fullName="fastqEnd2", shortName="o2", doc="Fastq end 2 out")
    public PrintStream o2;

    private Random rng;

    private Map<String, String> loadReference() {
        Map<String, String> ref = new TreeMap<String, String>();

        ReferenceSequence rseq;
        while ((rseq = REFERENCE.nextSequence()) != null) {
            String[] name = rseq.getName().split("\\s+");
            String seq = new String(rseq.getBases());

            ref.put(name[0], seq);
        }

        return ref;
    }

    private Map<String, int[]> loadReadStartDist(Map<String, String> ref) {
        Map<String, int[]> rs = new HashMap<String, int[]>();

        for (String chr : ref.keySet()) {
            int[] d = new int[ref.get(chr).length()];

            rs.put(chr, d);
        }

        TableReader tr = new TableReader(READ_START_DIST, new String[] { "chr", "start", "count" });

        for (Map<String, String> te : tr) {
            String chr = te.get("chr");
            int start = Integer.valueOf(te.get("start"));
            int count = Integer.valueOf(te.get("count"));

            rs.get(chr)[start] = count;
        }

        return rs;
    }

    private int[] loadDist(DataTable t, String binColumnName) {
        int maxBin = 0;

        for (Map<String, Object> o : t) {
            int bin = (Integer) o.get(binColumnName);
            //int count = (Integer) o.get("count");

            if (bin > maxBin) {
                maxBin = bin;
            }
        }

        int[] dist = new int[maxBin + 1];
        for (Map<String, Object> o : t) {
            int bin = (Integer) o.get(binColumnName);
            int count = (Integer) o.get("count");

            dist[bin] = count;
        }

        return dist;
    }

    private class FastEmpiricalDistribution {
        private double[] rates;
        private double[] cdf;
        private double tol = 1e-8;

        public FastEmpiricalDistribution(int[] dist) {
            this.rates = new double[dist.length];

            long sum = 0;
            for (int i = 0; i < dist.length; i++) {
                sum += dist[i];
            }

            for (int i = 0; i < dist.length; i++) {
                rates[i] = (double) dist[i] / (double) sum;
            }

            double cdfSum = 0.0;
            cdf = new double[rates.length];
            for (int i = 0; i < rates.length; i++) {
                cdfSum += rates[i];
                cdf[i] = cdfSum;
            }
        }

        private int findBin(double cdfValue) {
            int lowerBin = 0;
            int higherBin = cdf.length - 1;
            int middleBin = lowerBin + Math.round(((float) (higherBin - lowerBin))/2.0f);

            while (lowerBin != middleBin && middleBin != higherBin) {
                double lowerCdfValue = cdf[lowerBin];
                double middleCdfValue = cdf[middleBin];
                double higherCdfValue = cdf[higherBin];

                if (MoreMathUtils.equals(cdfValue, lowerCdfValue, tol))  {
                    return lowerBin;
                } else if (MoreMathUtils.equals(cdfValue, middleCdfValue, tol)) {
                    return middleBin;
                } else if (MoreMathUtils.equals(cdfValue, higherCdfValue, tol)) {
                    return higherBin;
                } else if (lowerCdfValue < cdfValue && cdfValue < middleCdfValue) {
                    higherBin = middleBin;
                    middleBin = lowerBin + Math.round(((float) (higherBin - lowerBin))/2.0f);
                } else if (middleCdfValue < cdfValue && cdfValue < higherCdfValue) {
                    lowerBin = middleBin;
                    middleBin = lowerBin + Math.round(((float) (higherBin - lowerBin))/2.0f);
                } else if (cdfValue < lowerCdfValue) {
                    return lowerBin;
                } else if (cdfValue > higherCdfValue) {
                    return higherBin;
                } else {
                    return middleBin;
                }
            }

            return middleBin;
        }

        public int draw() {
            return findBin(rng.nextDouble());
        }
    }

    private char getRandomNewBase(char oldBase) {
        final char[] bases = {'A', 'C', 'G', 'T'};

        char newBase = oldBase;
        while (newBase == oldBase) {
            int baseIndex = rng.nextInt(bases.length);
            newBase = bases[baseIndex];
        }

        return newBase;
    }

    private String generateRead(String fragment, int readLength, boolean isFirstEndOfRead, boolean fragmentIsNegativeStrand, FastEmpiricalDistribution errorNumRates, FastEmpiricalDistribution errorTypes, Map<String, FastEmpiricalDistribution> covRates) {
        StringBuilder read;
        if (isFirstEndOfRead) {
            read = new StringBuilder(fragment.substring(0, readLength + 10));
        } else {
            String fragmentRc = SequenceUtils.reverse(fragment);
            read = new StringBuilder(fragmentRc.substring(0, readLength + 10));
        }

        int numErrors = errorNumRates.draw();
        for (int e = 0; e < numErrors; e++) {
            int bin = errorTypes.draw();

            String type = "MM";
            switch (bin) {
                case 0: type = "MM"; break;
                case 1: type = "INS"; break;
                case 2: type = "DEL"; break;
            }

            String pkMarg = Joiner.on(".").join(type, isFirstEndOfRead, isFirstEndOfRead ? fragmentIsNegativeStrand : !fragmentIsNegativeStrand);
            int pos;
            do {
                pos = covRates.get(pkMarg).draw();
            } while (pos >= read.length());

            switch (bin) {
                case 0:
                    read.setCharAt(pos, getRandomNewBase(read.charAt(pos)));
                    break;
                case 1:
                    read.insert(pos, SequenceUtils.generateRandomNucleotideSequenceOfLengthN(rng.nextInt(2) + 1));
                    break;
                case 2:
                    int length = rng.nextInt(2) + 1;
                    read.delete(pos, pos + length);

                    break;
            }
        }

        return read.length() > readLength ? read.substring(0, readLength) : read.toString();
    }

    @Override
    public void execute() {
        if (SEED == null) { SEED = System.currentTimeMillis(); }
        rng = new Random(SEED);

        log.info("Loading data...");
        log.info("  reference sequence...");
        Map<String, String> ref = loadReference();

        log.info("  fragment start distribution...");
        Map<String, int[]> rs = loadReadStartDist(ref);

        log.info("  simulation profile...");
        DataTables ts = new DataTables(SIM_PROFILE);

        Map<String, Object> stats = ts.getTable("stats").iterator().next();

        int readLength = (Integer) stats.get("readLength");
        //int numReads = (Integer) stats.get("numReads");
        //int numPairs = (Integer) stats.get("numPairs");
        //int numChimericPairs = (Integer) stats.get("numChimericPairs");
        int numMismatches = (Integer) stats.get("numMismatches");
        int numInsertions = (Integer) stats.get("numInsertions");
        int numDeletions = (Integer) stats.get("numDeletions");

        log.info("    seed: {}", SEED);
        log.info("    readLength: {}", readLength);

        FastEmpiricalDistribution errorNumRates = new FastEmpiricalDistribution(loadDist(ts.getTable("errorNumDist"), "numErrorsPerRead"));
        FastEmpiricalDistribution fragmentSizeRates = new FastEmpiricalDistribution(loadDist(ts.getTable("fragmentSizeDist"), "fragmentSize"));
        FastEmpiricalDistribution errorTypes = new FastEmpiricalDistribution(new int[] { numMismatches, numInsertions, numDeletions });

        Map<String, int[]> covDists = new HashMap<String, int[]>();
        Map<String, FastEmpiricalDistribution> covRates = new HashMap<String, FastEmpiricalDistribution>();

        DataTable co = ts.getTable("covariateTable");
        for (Map<String, Object> ce : co) {
            String type = (String) ce.get("et");
            boolean isFirstEndOfRead = (Boolean) ce.get("isFirstEndOfRead");
            boolean isNegativeStrand = (Boolean) ce.get("isNegativeStrand");
            int position = (Integer) ce.get("position");
            int count = (Integer) ce.get("count");

            String pkMarg = Joiner.on(".").join(type, isFirstEndOfRead, isNegativeStrand);

            if (!covDists.containsKey(pkMarg)) { covDists.put(pkMarg, new int[readLength + 1]); }
            if (position >= 0 && position < readLength) { covDists.get(pkMarg)[position] = count; }
        }

        for (String pk : covDists.keySet()) {
            covRates.put(pk, new FastEmpiricalDistribution(covDists.get(pk)));
        }

        log.info("Simulating reads...");
        long fragmentIndex = 0;
        for (String chr : rs.keySet()) {
            log.info("  {}", chr);

            for (int i = 0; i < rs.get(chr).length; i++) {
                int numStarts = rs.get(chr)[i];

                for (int j = 0; j < numStarts; j++) {
                    int fragmentSize = 0;
                    while (fragmentSize <= readLength + 20) {
                        fragmentSize = fragmentSizeRates.draw();
                    }

                    if (i + fragmentSize < ref.get(chr).length()) {
                        String fragment = ref.get(chr).substring(i, i + fragmentSize);
                        boolean fragmentOnNegativeStrand = rng.nextBoolean();
                        if (fragmentOnNegativeStrand) {
                            fragment = SequenceUtils.reverseComplement(fragment);
                        }

                        String read1 = generateRead(fragment, readLength, true, fragmentOnNegativeStrand, errorNumRates, errorTypes, covRates);
                        String read2 = generateRead(fragment, readLength, false, fragmentOnNegativeStrand, errorNumRates, errorTypes, covRates);

                        String readName1 = String.format("@sim%d.%d %s:%d/1", SEED, fragmentIndex, chr, i);
                        String readName2 = String.format("@sim%d.%d %s:%d/2", SEED, fragmentIndex, chr, i);

                        o1.println(readName1);
                        o1.println(read1);
                        o1.println("+");
                        o1.println(StringUtils.repeat("I", read1.length()));

                        o2.println(readName2);
                        o2.println(read2);
                        o2.println("+");
                        o2.println(StringUtils.repeat("I", read2.length()));

                        fragmentIndex++;
                    }
                }
            }
        }
    }
}
