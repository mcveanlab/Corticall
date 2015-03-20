package uk.ac.ox.well.indiana.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTables;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.AlignmentUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

public class GenerateReadSimProfile extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public SAMFileReader BAM;

    @Output
    public PrintStream out;

    private enum ErrorType { MM, INS, DEL };

    private class ReadSimulationProfile {
        private Set<String> readNamesSeen = new HashSet<String>();
        private DataTables ts = new DataTables();

        private int readLength = 0;

        public ReadSimulationProfile() {
            ts.addTable("stats", "Basic stats from empirical data", "readLength", "numReads", "numReadsWithErrors", "numPairs", "numChimericPairs", "numMismatches", "numInsertions", "numDeletions");
            ts.addTable("fragmentSizeDist", "Distribution of fragment sizes", "fragmentSize", "count");
            ts.addTable("insertionSizeDist", "Distribution of insertion sizes", "insertionSize", "count");
            ts.addTable("deletionSizeDist", "Distribution of deletion sizes", "deletionSize", "count");
            ts.addTable("errorNumDist", "Distribution of number of errors in a read", "numErrors", "count");
            ts.addTable("covariateTable", "Table of error covariates", "type", "isFirstEndOfRead", "isNegativeStrand", "context", "position", "firstBaseOfError", "count");
        }

        public void accumulate(SAMRecord read) {
            ts.getTable("stats").increment("onlyrow", "numReads");

            if (readLength == 0 || read.getReadLength() > readLength) {
                readLength = read.getReadLength();

                ts.getTable("stats").set("onlyrow", "readLength", readLength);
            }

            if (!readNamesSeen.contains(read.getReadName())) {
                readNamesSeen.add(read.getReadName());
            } else {
                ts.getTable("stats").increment("onlyrow", "numPairs");

                if (!read.getReferenceName().equals(read.getMateReferenceName()) &&
                   (!read.getReadUnmappedFlag() || !read.getMateUnmappedFlag())) {
                    ts.getTable("stats").increment("onlyrow", "numChimericPairs");
                } else {
                    int fragmentSize = Math.abs(read.getInferredInsertSize());

                    if (fragmentSize > 0 && fragmentSize < 1000) {
                        String pk = String.format("%04d", fragmentSize);
                        ts.getTable("fragmentSizeDist").set(pk, "fragmentSize", fragmentSize);
                        ts.getTable("fragmentSizeDist").increment(pk, "count");
                    }
                }
            }

            Set<Integer> mismatchPositions = AlignmentUtils.getMismatchPositions(read);
            Set<Integer> insertionPositions = AlignmentUtils.getInsertionPositions(read);
            Set<Integer> deletionPositions = AlignmentUtils.getDeletionPositions(read);

            addErrors(read, mismatchPositions, ErrorType.MM);
            addErrors(read, insertionPositions, ErrorType.INS);
            addErrors(read, deletionPositions, ErrorType.DEL);

            ts.getTable("stats").add("onlyrow", "numMismatches", mismatchPositions.size());
            ts.getTable("stats").add("onlyrow", "numInsertions", insertionPositions.size());
            ts.getTable("stats").add("onlyrow", "numDeletions", deletionPositions.size());

            if (mismatchPositions.size() + insertionPositions.size() + deletionPositions.size() > 0) {
                ts.getTable("stats").increment("onlyrow", "numReadsWithErrors");
            }

            for (int size : AlignmentUtils.getInsertionSizes(read)) {
                ts.getTable("insertionSizeDist").set(String.format("%04d", size), "insertionSize", size);
                ts.getTable("insertionSizeDist").increment(String.format("%04d", size), "count");
            }

            for (int size : AlignmentUtils.getDeletionSizes(read)) {
                ts.getTable("deletionSizeDist").set(String.format("%04d", size), "deletionSize", size);
                ts.getTable("deletionSizeDist").increment(String.format("%04d", size), "count");
            }

            int numErrors = mismatchPositions.size() + insertionPositions.size() + deletionPositions.size();
            if (numErrors > 0) {
                ts.getTable("errorNumDist").set(String.format("%04d", numErrors), "numErrors", numErrors);
                ts.getTable("errorNumDist").increment(String.format("%04d", numErrors), "count");
            }
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

                if (!prevContext.contains("N") && pos >= 0 && pos < readLength && pos < readBases.length()) {
                    char firstBaseOfError = (type == ErrorType.DEL) ? '.' : read.getReadString().charAt(pos);
                    String pk = Joiner.on(".").join(type, read.getFirstOfPairFlag(), read.getReadNegativeStrandFlag(), prevContext, firstBaseOfError, String.format("%04d", pos));

                    ts.getTable("covariateTable").set(pk, "type", type);
                    ts.getTable("covariateTable").set(pk, "isFirstEndOfRead", read.getFirstOfPairFlag());
                    ts.getTable("covariateTable").set(pk, "isNegativeStrand", read.getReadNegativeStrandFlag());
                    ts.getTable("covariateTable").set(pk, "context", prevContext);
                    ts.getTable("covariateTable").set(pk, "firstBaseOfError", firstBaseOfError);
                    ts.getTable("covariateTable").set(pk, "position", pos);
                    ts.getTable("covariateTable").increment(pk, "count");
                }
            }
        }

        public void write(PrintStream out) {
            ts.write(out);
        }
    }

    @Override
    public void execute() {
        ReadSimulationProfile et = new ReadSimulationProfile();

        log.info("Processing reads...");
        String refName = "";
        for (SAMRecord read : BAM) {
            if (!read.getReferenceName().equals(refName)) {
                //if (refName.equals("Pf3D7_07_v3")) { break; }

                refName = read.getReferenceName();
                log.info("  {}", refName);
            }

            et.accumulate(read);
        }

        et.write(out);
    }
}
