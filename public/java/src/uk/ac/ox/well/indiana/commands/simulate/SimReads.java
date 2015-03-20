package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.ivy.util.StringUtils;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTable;
import uk.ac.ox.well.indiana.utils.containers.DataTables;
import uk.ac.ox.well.indiana.utils.containers.DataTree;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.statistics.distributions.EmpiricalDistribution;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class SimReads extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public FastaSequenceFile REFERENCE;

    @Argument(fullName="readStartDistribution", shortName="s", doc="Read start distribution")
    public File READ_START_DIST;

    @Argument(fullName="errorSimProfile", shortName="e", doc="Error simulation profile")
    public File SIM_PROFILE;

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

    private String generateRead(String fragment, int readLength, boolean isFirstEndOfRead, boolean fragmentIsNegativeStrand) {
        StringBuilder read;
        if (isFirstEndOfRead) {
            read = new StringBuilder(fragment);
        } else {
            String fragmentRc = SequenceUtils.reverseComplement(fragment);
            read = new StringBuilder(fragmentRc);
        }

        //boolean shouldHaveErrors = rateOfErrorfulReads.draw() == 1;
        boolean shouldHaveErrors = rng.nextDouble() <= rateOfErrorfulReads.getRates()[1];
        if (shouldHaveErrors) {
            int numErrors = errorNumRates.draw();

            int[] mmErrorRates = new int[readLength];
            int[] insErrorRates = new int[readLength];
            int[] delErrorRates = new int[readLength];

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
                insErrorRates[j - 2] = insCount;

                int delCount = 0;
                if (covTree.has("DEL", isFirstEndOfRead, fragmentIsNegativeStrand, context, j - 2, ".")) {
                    delCount += (Integer) covTree.get("DEL", isFirstEndOfRead, fragmentIsNegativeStrand, context, j - 2, ".").iterator().next();
                }
                delErrorRates[j - 2] = delCount;

                bases.remove(currentBase);

                int mmCount = 0;
                for (String base : bases) {
                    if (covTree.has("MM", isFirstEndOfRead, fragmentIsNegativeStrand, context, j - 2, base)) {
                        mmCount += (Integer) covTree.get("MM", isFirstEndOfRead, fragmentIsNegativeStrand, context, j - 2, base).iterator().next();
                    }
                }
                mmErrorRates[j - 2] = mmCount;
            }

            EmpiricalDistribution mmd  = new EmpiricalDistribution(mmErrorRates);
            EmpiricalDistribution insd = new EmpiricalDistribution(insErrorRates);
            EmpiricalDistribution deld = new EmpiricalDistribution(delErrorRates);

            String[] bases = new String[] { "A", "C", "G", "T" };

            for (int i = 0; i < numErrors; i++) {
                int errorTypeBin = errorTypes.draw();

                if (errorTypeBin == 0) { // mismatches
                    int errorPos = mmd.draw() + 2;
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

                    EmpiricalDistribution bd = new EmpiricalDistribution(baseDist);
                    String newBase = bases[bd.draw()];

                    read.setCharAt(errorPos, newBase.charAt(0));
                } else if (errorTypeBin == 1) { // insertions
                    int errorPos = insd.draw() + 2;
                    String context = read.substring(errorPos - 2, errorPos);

                    int[] baseDist = new int[4];
                    for (int c = 0; c < bases.length; c++) {
                        if (covTree.has("INS", isFirstEndOfRead, fragmentIsNegativeStrand, context, errorPos - 2, bases[c])) {
                            baseDist[c] = (Integer) covTree.get("INS", isFirstEndOfRead, fragmentIsNegativeStrand, context, errorPos - 2, bases[c]).iterator().next();
                        }
                    }

                    EmpiricalDistribution bd = new EmpiricalDistribution(baseDist);
                    String newBase = bases[bd.draw()];

                    int insertionLength = insertionSizeRates.draw() - 1;
                    String restOfInsertion = new String(SequenceUtils.generateRandomNucleotideSequenceOfLengthN(insertionLength));

                    String insertionSeq = newBase + restOfInsertion;
                    read.insert(errorPos, insertionSeq);
                } else if (errorTypeBin == 2) { // deletion
                    int errorPos = deld.draw() + 2;
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

        fragmentSizeRates = new EmpiricalDistribution(loadDist(ts.getTable("fragmentSizeDist"), "fragmentSize"));
        rateOfChimeras = new EmpiricalDistribution(new int[] { numPairs - numChimericPairs, numChimericPairs });
        rateOfErrorfulReads = new EmpiricalDistribution(new int[] { numReads - numReadsWithErrors, numReadsWithErrors });
        errorNumRates = new EmpiricalDistribution(loadDist(ts.getTable("errorNumDist"), "numErrors"));
        errorTypes = new EmpiricalDistribution(new int[] { numMismatches, numInsertions, numDeletions });
        insertionSizeRates = new EmpiricalDistribution(loadDist(ts.getTable("insertionSizeDist"), "insertionSize"));
        deletionSizeRates = new EmpiricalDistribution(loadDist(ts.getTable("deletionSizeDist"), "deletionSize"));

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

    @Override
    public void execute() {
        if (SEED == null) { SEED = System.currentTimeMillis(); }
        rng = new Random(SEED);

        log.info("Loading data...");
        log.info("  simulation profile...");
        DataTables ts = new DataTables(SIM_PROFILE);
        initializeEmpiricalDistributions(ts);
        int readLength = (Integer) ts.getTable("stats").iterator().next().get("readLength");

        log.info("    seed: {}", SEED);
        log.info("    readLength: {}", readLength);

        log.info("  reference sequence...");
        Map<String, String> ref = loadReference();
        chrNames.addAll(ref.keySet());

        log.info("  fragment start distribution...");
        Map<String, int[]> rs = loadReadStartDist(ref);

        log.info("Simulating reads...");
        long fragmentIndex = 0;
        for (String chr : ref.keySet()) {
            if (chr.equals(CHR)) {
                log.info("  {}", chr);

                for (int i = 0; i < rs.get(chr).length; i++) {
                    if (i % (rs.get(chr).length / 10) == 0) {
                        log.info("    {}/{} bp", i, rs.get(chr).length);
                    }

                    int numStarts = rs.get(chr)[i];

                    for (int j = 0; j < numStarts; j++) {
                        String fragment = generateFragment(ref, chr, i, readLength);

                        if (!fragment.isEmpty()) {
                            boolean fragmentOnNegativeStrand = rng.nextBoolean();
                            if (fragmentOnNegativeStrand) {
                                fragment = SequenceUtils.reverseComplement(fragment);
                            }

                            String read1 = generateRead(fragment, readLength, true, fragmentOnNegativeStrand);
                            String read2 = generateRead(fragment, readLength, false, fragmentOnNegativeStrand);

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
    }
}
