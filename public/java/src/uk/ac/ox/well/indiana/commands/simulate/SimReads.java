package uk.ac.ox.well.indiana.commands.simulate;

import cern.jet.random.Empirical;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import org.apache.ivy.util.StringUtils;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTable;
import uk.ac.ox.well.indiana.utils.containers.DataTables;
import uk.ac.ox.well.indiana.utils.containers.DataTree;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.statistics.distributions.EmpiricalDistribution;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class SimReads extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public FastaSequenceFile REFERENCE;

    @Argument(fullName="readProfile", shortName="s", doc="Read simulation profile")
    public File READ_PROFILE;

    @Argument(fullName="errorSimProfile", shortName="e", doc="Error simulation profile")
    public File ERROR_PROFILE;

    @Argument(fullName="repetitiveRegions", shortName="b", doc="Repetitive regions (bed)")
    public File REPETITIVE_REGIONS;

    @Argument(fullName="chr", shortName="c", doc="Chr to process")
    public String CHR;

    @Argument(fullName="seed", shortName="g", doc="Seed for RNG", required=false)
    public Long SEED;

    @Output(fullName="fastqEnd1", shortName="o1", doc="q end 1 out")
    public PrintStream o1;

    @Output(fullName="fastqEnd2", shortName="o2", doc="q end 2 out")
    public PrintStream o2;

    private Random rng;
    private EmpiricalDistribution fragmentSizeRates;
    private EmpiricalDistribution rateOfChimeras;
    private EmpiricalDistribution rateOfErrorfulReads;
    private EmpiricalDistribution errorNumRates;
    private EmpiricalDistribution errorTypes;
    private EmpiricalDistribution insertionSizeRates;
    private EmpiricalDistribution deletionSizeRates;
    private DataTree covTree;

    private List<String> chrNames = new ArrayList<String>();
    private int numChimerasIntroduced = 0;

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

    private class Entry {
        public int numReadsWithErrors;
        public int numReads;
        public int numFragmentsWithErrors;
        public int numFragments;
    }

    private Map<String, Map<Integer, Entry>> loadReadProfile() {
        Map<String, Map<Integer, Entry>> table = new HashMap<String, Map<Integer, Entry>>();

        LineReader lr = new LineReader(READ_PROFILE);
        lr.getNextRecord();

        while (lr.hasNext()) {
            String line = lr.getNextRecord();
            String[] fields = line.split("\\s+");

            String chr = fields[0];
            int pos = Integer.valueOf(fields[1]);

            Entry e = new Entry();
            e.numReadsWithErrors = Integer.valueOf(fields[2]);
            e.numReads = Integer.valueOf(fields[3]);
            e.numFragmentsWithErrors = Integer.valueOf(fields[4]);
            e.numFragments = Integer.valueOf(fields[5]);

            if (!table.containsKey(chr)) {
                table.put(chr, new HashMap<Integer, Entry>());
            }

            if (!table.get(chr).containsKey(pos)) {
                table.get(chr).put(pos, e);
            }
        }

        return table;
    }

    private int[] loadDist(DataTable t, String binColumnName) {
        int maxBin = 0;

        for (Map<String, Object> o : t) {
            if (!o.containsKey(binColumnName)) {
                throw new IndianaException("Column " + binColumnName + " does not exist in DataTable " + t.getTableName());
            }

            int bin = (Integer) o.get(binColumnName);

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

    private Map<String, IntervalTreeMap<Integer>> loadRepetitiveRegions() {
        Map<String, IntervalTreeMap<Integer>> rr = new HashMap<String, IntervalTreeMap<Integer>>();
        TableReader tr = new TableReader(REPETITIVE_REGIONS, new String[] { "sequence", "begin", "end" });

        for (Map<String, String> te : tr) {
            String chr = te.get("sequence");
            int begin = Integer.valueOf(te.get("begin"));
            int end = Integer.valueOf(te.get("end"));

            Interval interval = new Interval(chr, begin, end);

            if (!rr.containsKey(chr)) {
                rr.put(chr, new IntervalTreeMap<Integer>());
            }

            rr.get(chr).put(interval, null);
        }

        return rr;
    }

    private int computeNumErrors(double[] posErrorRates, int readLength) {
        int numErrors = 0;
        for (int i = 2; i < readLength + 2; i++) {
            double val = rng.nextDouble();

            if (posErrorRates.length > i && val <= posErrorRates[i]) {
                numErrors++;
            }
        }

        return numErrors;
    }

    private String generateRead(String fragment, double[] posErrorRates, int readLength, boolean isFirstEndOfRead, boolean fragmentIsNegativeStrand, boolean shouldHaveErrors) {
        StringBuilder read;
        if (isFirstEndOfRead) {
            read = new StringBuilder(fragment);
        } else {
            String fragmentRc = SequenceUtils.reverseComplement(fragment);
            read = new StringBuilder(fragmentRc);

            double[] posErrorRatesRc = new double[posErrorRates.length];
            for (int i = 0; i < posErrorRates.length; i++) {
                posErrorRatesRc[i] = posErrorRates[posErrorRates.length - i - 1];
            }

            posErrorRates = posErrorRatesRc;
        }

        double[] posErrorNorm = new double[readLength];
        double posErrorSum = 0;
        for (int i = 0; i < readLength; i++) {
            if (i + 2 < posErrorRates.length) {
                posErrorSum += posErrorRates[i + 2];
            }
        }
        for (int i = 0; i < readLength; i++) {
            if (i + 2 < posErrorRates.length) {
                posErrorNorm[i] = posErrorRates[i + 2] / posErrorSum;
            }
        }

        EmpiricalDistribution posd = new EmpiricalDistribution(posErrorRates, rng);
        //EmpiricalDistribution psnd = new EmpiricalDistribution(posErrorNorm, rng);

        if (shouldHaveErrors) {
            int numErrors = computeNumErrors(posErrorRates, readLength);

            double[] insErrorDist = new double[readLength];
            double insErrorSum = 0;

            double[] delErrorDist = new double[readLength];
            double delErrorSum = 0;

            double[] mmErrorDist = new double[readLength];
            double mmErrorSum = 0;

            for (int j = 2; j < readLength + 2; j++) {
                String context = read.substring(j - 2, j);

                String currentBase = String.valueOf(read.charAt(j));
                Set<String> bases = new HashSet<String>(Arrays.asList("A", "C", "G", "T"));

                int insCount = 0;
                for (String base : bases) {
                    if (covTree.has("INS", isFirstEndOfRead, fragmentIsNegativeStrand, context, j - 2, base)) {
                        insCount += (Integer) covTree.get("INS", isFirstEndOfRead, fragmentIsNegativeStrand, context, j - 2, base).iterator().next();
                    }
                }
                insErrorDist[j - 2] = insCount;
                insErrorSum += insCount;

                int delCount = 0;
                if (covTree.has("DEL", isFirstEndOfRead, fragmentIsNegativeStrand, context, j - 2, ".")) {
                    delCount += (Integer) covTree.get("DEL", isFirstEndOfRead, fragmentIsNegativeStrand, context, j - 2, ".").iterator().next();
                }
                delErrorDist[j - 2] = delCount;
                delErrorSum += delCount;

                bases.remove(currentBase);

                int mmCount = 0;
                for (String base : bases) {
                    if (covTree.has("MM", isFirstEndOfRead, fragmentIsNegativeStrand, context, j - 2, base)) {
                        mmCount += (Integer) covTree.get("MM", isFirstEndOfRead, fragmentIsNegativeStrand, context, j - 2, base).iterator().next();
                    }
                }
                mmErrorDist[j - 2] = mmCount;
                mmErrorSum += mmCount;
            }

            double insErrorNorm = 0;
            double delErrorNorm = 0;
            double mmErrorNorm = 0;
            for (int i = 0; i < readLength; i++) {
                insErrorDist[i] = (posErrorNorm[i]) * (insErrorDist[i] / insErrorSum);
                delErrorDist[i] = (posErrorNorm[i]) * (delErrorDist[i] / delErrorSum);
                mmErrorDist[i]  = (posErrorNorm[i]) * (mmErrorDist[i]  / mmErrorSum);

                insErrorNorm += insErrorDist[i];
                delErrorNorm += delErrorDist[i];
                mmErrorNorm += mmErrorDist[i];
            }

            for (int i = 0; i < readLength; i++) {
                insErrorDist[i] /= insErrorNorm;
                delErrorDist[i] /= delErrorNorm;
                mmErrorDist[i]  /= mmErrorNorm;
            }

            EmpiricalDistribution mmd  = new EmpiricalDistribution(mmErrorDist, rng);
            EmpiricalDistribution insd = new EmpiricalDistribution(insErrorDist, rng);
            EmpiricalDistribution deld = new EmpiricalDistribution(delErrorDist, rng);

            String[] bases = new String[] { "A", "C", "G", "T" };

            for (int i = 0; i < numErrors; i++) {
                int errorTypeBin = errorTypes.draw();

                if (errorTypeBin == 0) { // mismatches
                    int errorPos;
                    do {
                        errorPos = mmd.draw() + 2;
                    } while (errorPos > read.length());
                    String correctBase = String.valueOf(read.charAt(errorPos));
                    String context = read.substring(errorPos - 2, errorPos);

                    int[] baseDist = new int[4];
                    for (int c = 0; c < bases.length; c++) {
                        if (!bases[c].equals(correctBase)) {
                            if (covTree.has("MM", isFirstEndOfRead, fragmentIsNegativeStrand, context, errorPos - 2, bases[c])) {
                                baseDist[c] = (Integer) covTree.get("MM", isFirstEndOfRead, fragmentIsNegativeStrand, context, errorPos - 2, bases[c]).iterator().next();
                            }
                        }
                    }

                    EmpiricalDistribution bd = new EmpiricalDistribution(baseDist, rng);
                    String newBase = bases[bd.draw()];

                    read.setCharAt(errorPos, newBase.charAt(0));
                } else if (errorTypeBin == 1) { // insertions
                    int errorPos;
                    do {
                        errorPos = insd.draw() + 2;
                    } while (errorPos > read.length());
                    String context = read.substring(errorPos - 2, errorPos);

                    int[] baseDist = new int[4];
                    for (int c = 0; c < bases.length; c++) {
                        if (covTree.has("INS", isFirstEndOfRead, fragmentIsNegativeStrand, context, errorPos - 2, bases[c])) {
                            baseDist[c] = (Integer) covTree.get("INS", isFirstEndOfRead, fragmentIsNegativeStrand, context, errorPos - 2, bases[c]).iterator().next();
                        }
                    }

                    EmpiricalDistribution bd = new EmpiricalDistribution(baseDist, rng);
                    String newBase = bases[bd.draw()];

                    int insertionLength = insertionSizeRates.draw() - 1;
                    String restOfInsertion = new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(insertionLength));

                    String insertionSeq = newBase + restOfInsertion;
                    read.insert(errorPos, insertionSeq);
                } else if (errorTypeBin == 2) { // deletion
                    int errorPos;
                    do {
                        errorPos = deld.draw() + 2;
                    } while (errorPos > read.length());
                    int deletionLength = deletionSizeRates.draw();

                    read.delete(errorPos, errorPos + deletionLength);
                }
            }

        }

        return (read.length() > readLength + 2) ? read.substring(2, readLength + 2) : null;
    }

    private void initializeEmpiricalDistributions(DataTables ts) {
        Map<String, Object> stats = ts.getTable("stats").iterator().next();

        int numReads = (Integer) stats.get("numReads");
        int numReadsWithErrors = (Integer) stats.get("numReadsWithErrors");
        int numPairs = (Integer) stats.get("numPairs");
        int numChimericPairs = (Integer) stats.get("numChimericPairs");
        int numMismatches = (Integer) stats.get("numMismatches");
        int numInsertions = (Integer) stats.get("numInsertions");
        int numDeletions = (Integer) stats.get("numDeletions");

        fragmentSizeRates = new EmpiricalDistribution(loadDist(ts.getTable("fragmentSizeDist"), "fragmentSize"), rng);
        rateOfChimeras = new EmpiricalDistribution(new int[] { numPairs - numChimericPairs, numChimericPairs }, rng);
        rateOfErrorfulReads = new EmpiricalDistribution(new int[] { numReads - numReadsWithErrors, numReadsWithErrors }, rng);
        errorNumRates = new EmpiricalDistribution(loadDist(ts.getTable("errorNumDist"), "numErrors"), rng);
        errorTypes = new EmpiricalDistribution(new int[] { numMismatches, numInsertions, numDeletions }, rng);
        insertionSizeRates = new EmpiricalDistribution(loadDist(ts.getTable("insertionSizeDist"), "insertionSize"), rng);
        deletionSizeRates = new EmpiricalDistribution(loadDist(ts.getTable("deletionSizeDist"), "deletionSize"), rng);

        covTree = new DataTree();
        for (Map<String, Object> e : ts.getTable("covariateTable")) {
            String type = (String) e.get("type");
            boolean isFirstEndOfRead = (Boolean) e.get("isFirstEndOfRead");
            boolean isNegativeStrand = (Boolean) e.get("isNegativeStrand");
            String context = (String) e.get("context");
            int position = (Integer) e.get("position");
            String firstBaseOfError = (String) e.get("firstBaseOfError");
            int count = (Integer) e.get("count");

            covTree.add(type, isFirstEndOfRead, isNegativeStrand, context, position, firstBaseOfError, count);
        }
    }

    private String generateFragment(Map<String, String> ref, String chr, int pos, int readLength) {
        int fragmentSize = 0;
        while (fragmentSize <= readLength + 4) {
            fragmentSize = fragmentSizeRates.draw();
        }

        String fragment = "";
        if (pos + fragmentSize < ref.get(chr).length()) {
            fragment = ref.get(chr).substring(pos, pos + fragmentSize);

            boolean shouldBeChimeric = rng.nextDouble() <= rateOfChimeras.getRates()[1];
            if (shouldBeChimeric) {
                numChimerasIntroduced++;
                String chimericChr = chr;

                while (chimericChr.equals(chr)) {
                    int chrIndex = rng.nextInt(chrNames.size());
                    chimericChr = chrNames.get(chrIndex);
                }

                int chimericPos = rng.nextInt(ref.get(chimericChr).length() - 3*readLength);
                String chimericFragment = ref.get(chimericChr).substring(chimericPos, chimericPos + 2*readLength);

                fragment += chimericFragment;
            }
        }

        return fragment;
    }

    private double[] getPositionalErrorRates(Map<Integer, Entry> es, int pos, int length) {
        double[] er = new double[length];

        for (int i = 0; i < length; i++) {
            if (es.containsKey(pos + i)) {
                er[i] = (double) es.get(pos + i).numReadsWithErrors / (double) es.get(pos + i).numReads;
            }
        }

        return er;
    }

    @Override
    public void execute() {
        if (SEED == null) { SEED = System.currentTimeMillis(); }
        rng = new Random(SEED);

        log.info("Loading data...");
        log.info("  error profile...");
        DataTables ts = new DataTables(ERROR_PROFILE);
        initializeEmpiricalDistributions(ts);
        int readLength = (Integer) ts.getTable("stats").iterator().next().get("readLength");

        log.info("    seed: {}", SEED);
        log.info("    readLength: {}", readLength);

        log.info("  read profile...");
        Map<String, Map<Integer, Entry>> rs = loadReadProfile();

        log.info("  reference sequence...");
        Map<String, String> ref = loadReference();
        chrNames.addAll(ref.keySet());

        log.info("  repetitive regions...");
        Map<String, IntervalTreeMap<Integer>> rr = loadRepetitiveRegions();

        log.info("Simulating reads...");
        long fragmentIndex = 0;
        int alteredFragments = 0;
        for (String chr : rs.keySet()) {
            if (chr.equals(CHR) && ref.containsKey(chr)) {
                log.info("  {}", chr);

                for (int i = 0; i < ref.get(chr).length(); i++) {
                    if (i % (ref.get(chr).length() / 10) == 0) {
                        log.info("    {}/{} bp", i, ref.get(chr).length());
                    }

                    int numStarts = rs.get(chr).containsKey(i) ? rs.get(chr).get(i).numFragments : 0;
                    int numErrors = rs.get(chr).containsKey(i) ? rs.get(chr).get(i).numFragmentsWithErrors : 0;

                    for (int j = 0; j < numStarts; j++) {
                        String fragment = generateFragment(ref, chr, i, readLength);
                        double[] posErrorRates = getPositionalErrorRates(rs.get(chr), i, fragment.length());

                        if (!fragment.isEmpty()) {
                            boolean fragmentOnNegativeStrand = rng.nextBoolean();

                            Interval interval = new Interval(chr, i, i + fragment.length());
                            if (rr.get(chr).containsContained(interval) && rng.nextDouble() <= 0.50) {
                                int deletionPos = rng.nextInt(fragment.length());
                                int deletionLength = 5 + rng.nextInt(20);

                                if (deletionPos + deletionLength < fragment.length() && fragment.length() - deletionLength > readLength + 20) {
                                    //log.info("{} {} {}", deletionPos, deletionLength, fragment.length());

                                    StringBuilder sb = new StringBuilder(fragment);
                                    sb.delete(deletionPos, deletionPos + deletionLength);

                                    fragment = sb.toString();

                                    alteredFragments++;
                                }
                            }

                            if (fragmentOnNegativeStrand) {
                                fragment = SequenceUtils.reverseComplement(fragment);
                            }


                            String read1 = generateRead(fragment, posErrorRates, readLength, true, fragmentOnNegativeStrand, numErrors < j);
                            String read2 = generateRead(fragment, posErrorRates, readLength, false, fragmentOnNegativeStrand, numErrors < j);

                            if (read1 != null && read2 != null) {
                                String readName1 = String.format("@sim%d.%s.%d %s:%d/1", SEED, chr, fragmentIndex, chr, i);
                                String readName2 = String.format("@sim%d.%s.%d %s:%d/2", SEED, chr, fragmentIndex, chr, i);

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

        log.info("Stats");
        log.info("  paired-end reads: {}", fragmentIndex);
        log.info("  chimeric pairs:   {}", numChimerasIntroduced);
        log.info("  altered fragment: {}", alteredFragments);
    }
}
