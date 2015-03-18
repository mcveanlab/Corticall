package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GenerateReadSimProfile extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public SAMFileReader BAM;

    @Output
    public PrintStream out;

    private enum ErrorType { MM, INS, DEL };

    private class CovariateTable {
        private Map<ErrorType, Map<Boolean, Map<Boolean, Map<String, Map<Integer, Long>>>>> table = new TreeMap<ErrorType, Map<Boolean, Map<Boolean, Map<String, Map<Integer, Long>>>>>();

        public void increment(ErrorType type, boolean isFirstEndOfRead, boolean isNegativeStrand, String context, int position) {
            if (!table.containsKey(type)) {
                table.put(type, new TreeMap<Boolean, Map<Boolean, Map<String, Map<Integer, Long>>>>());
            }

            if (!table.get(type).containsKey(isFirstEndOfRead)) {
                table.get(type).put(isFirstEndOfRead, new TreeMap<Boolean, Map<String, Map<Integer, Long>>>());
            }

            if (!table.get(type).get(isFirstEndOfRead).containsKey(isNegativeStrand)) {
                table.get(type).get(isFirstEndOfRead).put(isNegativeStrand, new TreeMap<String, Map<Integer, Long>>());
            }

            if (!table.get(type).get(isFirstEndOfRead).get(isNegativeStrand).containsKey(context)) {
                table.get(type).get(isFirstEndOfRead).get(isNegativeStrand).put(context, new TreeMap<Integer, Long>());
            }

            if (!table.get(type).get(isFirstEndOfRead).get(isNegativeStrand).get(context).containsKey(position)) {
                table.get(type).get(isFirstEndOfRead).get(isNegativeStrand).get(context).put(position, 0l);
            }

            long oldValue = table.get(type).get(isFirstEndOfRead).get(isNegativeStrand).get(context).get(position);
            table.get(type).get(isFirstEndOfRead).get(isNegativeStrand).get(context).put(position, oldValue + 1);
        }

        public String toString() {
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            PrintStream ps = new PrintStream(baos);

            TableWriter tw = new TableWriter(ps);

            for (ErrorType et : table.keySet()) {
                for (Boolean isFirstEndOfRead : table.get(et).keySet()) {
                    for (Boolean isNegativeStrand : table.get(et).get(isFirstEndOfRead).keySet()) {
                        for (String context : table.get(et).get(isFirstEndOfRead).get(isNegativeStrand).keySet()) {
                            for (Integer position : table.get(et).get(isFirstEndOfRead).get(isNegativeStrand).get(context).keySet()) {
                                long count = table.get(et).get(isFirstEndOfRead).get(isNegativeStrand).get(context).get(position);

                                Map<String, String> te = new LinkedHashMap<String, String>();

                                te.put("covariateTable", "covariateTable");
                                te.put("et", et.name());
                                te.put("isFirstEndOfRead", String.valueOf(isFirstEndOfRead));
                                te.put("isNegativeStrand", String.valueOf(isNegativeStrand));
                                te.put("context", context);
                                te.put("position", String.valueOf(position));
                                te.put("count", String.valueOf(count));

                                tw.addEntry(te);
                            }
                        }
                    }
                }
            }

            return baos.toString();
        }
    }

    private class Histogram {
        private Map<Integer, Integer> hist = new TreeMap<Integer, Integer>();
        private String tableName;
        private String keyName;
        private String valName;

        public Histogram(String tableName, String keyName, String valName) {
            this.tableName = tableName;
            this.keyName = keyName;
            this.valName = valName;
        }

        public void increment(int bin) {
            if (!hist.containsKey(bin)) { hist.put(bin, 0); }

            hist.put(bin, hist.get(bin) + 1);
        }

        public String toString() {
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            PrintStream ps = new PrintStream(baos);

            TableWriter tw = new TableWriter(ps);

            for (int bin : hist.keySet()) {
                Map<String, String> te = new LinkedHashMap<String, String>();
                te.put(tableName, tableName);
                te.put(keyName, String.valueOf(bin));
                te.put(valName, String.valueOf(hist.get(bin)));

                tw.addEntry(te);
            }

            return baos.toString();
        }
    }

    private class ReadSimulationProfile {
        private int readLength = 0;
        private long numReads = 0;
        private long numPairs = 0;
        private long numChimericPairs = 0;

        private int numInsertions = 0;
        private int numDeletions = 0;
        private int numMismatches = 0;

        private Set<String> readNamesSeen = new HashSet<String>();

        private CovariateTable t = new CovariateTable();
        private Histogram fs = new Histogram("fragmentSizeDist", "fragmentSize", "count");
        private Histogram es = new Histogram("errorNumDist", "numErrorsPerRead", "count");

        private Pattern p = Pattern.compile("[0-9]+|[A-Z]|\\^[A-Z]+");

        private int computeInitialPosition(SAMRecord read) {
            int currentPosition = 0;

            if (read.getCigar().getCigarElements().size() > 0 && read.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.S)) {
                currentPosition = read.getCigar().getCigarElement(0).getLength();
            }

            return currentPosition;
        }

        public void accumulate(SAMRecord read) {
            numReads++;

            if (readLength == 0 || read.getReadLength() > readLength) {
                readLength = read.getReadLength();
            }

            if (!readNamesSeen.contains(read.getReadName())) {
                readNamesSeen.add(read.getReadName());
            } else {
                numPairs++;
                if (!read.getReferenceName().equals(read.getMateReferenceName()) &&
                   (!read.getReadUnmappedFlag() || !read.getMateUnmappedFlag())) {
                    numChimericPairs++;
                } else {
                    int fragmentSize = Math.abs(read.getInferredInsertSize());

                    if (fragmentSize > 0 && fragmentSize < 1000) {
                        fs.increment(fragmentSize);
                    }
                }
            }

            String md = read.getStringAttribute("MD");

            int initialPosition = computeInitialPosition(read);
            int currentPosition = initialPosition;

            Set<Integer> insertionPositions = new TreeSet<Integer>();
            Set<Integer> deletionPositions = new TreeSet<Integer>();
            Set<Integer> mismatchPositions = new TreeSet<Integer>();
            Set<Integer> errorPositions = new TreeSet<Integer>();

            for (CigarElement ce : read.getCigar().getCigarElements()) {
                if (!ce.getOperator().equals(CigarOperator.H) && !ce.getOperator().equals(CigarOperator.S)) {
                    if (ce.getOperator().equals(CigarOperator.I)) {
                        insertionPositions.add(currentPosition);
                        errorPositions.add(currentPosition);

                        numInsertions++;
                    } else if (ce.getOperator().equals(CigarOperator.D)) {
                        deletionPositions.add(currentPosition);
                        errorPositions.add(currentPosition);

                        numDeletions++;
                    }

                    currentPosition += ce.getLength();
                }
            }

            addErrors(read, insertionPositions, ErrorType.INS);
            addErrors(read, deletionPositions, ErrorType.DEL);

            if (md != null) {
                Matcher m = p.matcher(md);

                currentPosition = initialPosition;
                if (read.getCigar().getCigarElements().size() > 0 && read.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.S)) {
                    currentPosition = read.getCigar().getCigarElement(0).getLength();
                }

                while (m.find()) {
                    String piece = md.substring(m.start(), m.end());

                    if (!piece.contains("^")) {
                        if (piece.matches("[A-Z]")) {
                            numMismatches++;
                            mismatchPositions.add(currentPosition);
                            errorPositions.add(currentPosition);

                            currentPosition++;
                        } else {
                            currentPosition += Integer.parseInt(piece);
                        }
                    }
                }

                addErrors(read, mismatchPositions, ErrorType.MM);
            }

            es.increment(errorPositions.size());
        }

        private Set<Integer> flipStrand(Set<Integer> positions, int readLength) {
            Set<Integer> flippedPositions = new TreeSet<Integer>();

            for (int pos : positions) {
                flippedPositions.add(readLength - pos - 1);
            }

            return flippedPositions;
        }

        public void addErrors(SAMRecord read, Set<Integer> positions, ErrorType type) {
            String readBases = read.getReadString();
            if (read.getReadNegativeStrandFlag()) {
                readBases = SequenceUtils.reverseComplement(readBases);
                positions = flipStrand(positions, read.getReadLength());
            }

            for (int pos : positions) {
                String prevContext = "NN";

                if (pos <= 1 && read.getAlignmentStart() > 3) {
                    int start = read.getAlignmentStart() + pos;

                    prevContext = new String(REFERENCE.getSubsequenceAt(read.getReferenceName(), start - 2, start - 1).getBases());
                } else if (pos >= read.getReadLength() - 2 && read.getAlignmentEnd() < REFERENCE.getSequence(read.getReferenceName()).length() - 3) {
                    int start = read.getAlignmentEnd() - 1;

                    prevContext = new String(REFERENCE.getSubsequenceAt(read.getReferenceName(), start - 1, start).getBases());
                } else if (pos >= 2 && pos <= read.getReadLength() - 3) {
                    prevContext = readBases.substring(pos - 2, pos);
                }

                if (!prevContext.contains("N")) {
                    t.increment(type, read.getFirstOfPairFlag(), read.getReadNegativeStrandFlag(), prevContext, pos);
                }
            }
        }

        public String toString() {
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            PrintStream ps = new PrintStream(baos);

            TableWriter tw = new TableWriter(ps);

            Map<String, String> te = new LinkedHashMap<String, String>();
            te.put("stats", "stats");
            te.put("readLength", String.valueOf(readLength));
            te.put("numReads", String.valueOf(numReads));
            te.put("numPairs", String.valueOf(numPairs));
            te.put("numChimericPairs", String.valueOf(numChimericPairs));
            te.put("numMismatches", String.valueOf(numMismatches));
            te.put("numInsertions", String.valueOf(numInsertions));
            te.put("numDeletions", String.valueOf(numDeletions));

            tw.addEntry(te);

            return baos.toString() + "\n" +
                   es.toString()   + "\n" +
                   fs.toString()   + "\n" +
                   t.toString();
        }
    }

    @Override
    public void execute() {
        ReadSimulationProfile et = new ReadSimulationProfile();

        log.info("Processing reads...");
        String refName = "";
        for (SAMRecord read : BAM) {
            if (!read.getReferenceName().equals(refName)) {
                //if (refName.equals("Pf3D7_03_v3")) { break; }

                refName = read.getReferenceName();
                log.info("  {}", refName);
            }

            et.accumulate(read);
        }

        out.println(et);
    }
}
