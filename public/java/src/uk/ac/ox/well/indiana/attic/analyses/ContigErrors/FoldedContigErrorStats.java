package uk.ac.ox.well.indiana.attic.analyses.ContigErrors;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.PrintStream;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class FoldedContigErrorStats extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (in BAM format)")
    public SAMFileReader CONTIGS;

    @Output
    public PrintStream out;

    @Output(fullName="out2", shortName="o2")
    public PrintStream out2;

    private void incrementErrorsArray(Set<Integer> errorSet, int[] errors, int[] norms, int contigLength) {
        for (Integer error : errorSet) {
            if (error < contigLength/2) {
                errors[error]++;
            } else {
                if (contigLength - error >= 0) {
                    errors[contigLength - error]++;
                }
            }
        }

        for (int i = 0; i < contigLength; i++) {
            if (i < contigLength/2) {
                norms[i]++;
            } else {
                if (contigLength - i >= 0) {
                    norms[contigLength - i]++;
                }
            }
        }
    }

    @Override
    public void execute() {
        Pattern p = Pattern.compile("[0-9]+|[A-Z]|\\^[A-Z]+");

        TableWriter tw2 = new TableWriter(out2);

        int lengthMax = 100000;

        int[] mismatches = new int[lengthMax];
        int[] mismatchesTotal = new int[lengthMax];
        int[] insertions = new int[lengthMax];
        int[] insertionsTotal = new int[lengthMax];
        int[] deletions = new int[lengthMax];
        int[] deletionsTotal = new int[lengthMax];
        int[] all = new int[lengthMax];
        int[] allTotal = new int[lengthMax];

        for (SAMRecord read : CONTIGS) {
            if (read.getAlignmentStart() > 0) { // && read.getReadLength() >= 0 && read.getReadLength() < lengthMax) {
                String md = read.getStringAttribute("MD");

                int numMatches = 0;
                int numMismatches = 0;
                int numDeletions = 0;
                int numInsertions = 0;

                int currentPosition = 0;
                Set<Integer> insertionPositions = new TreeSet<Integer>();
                Set<Integer> deletionPositions = new TreeSet<Integer>();
                Set<Integer> mismatchPositions = new TreeSet<Integer>();
                Set<Integer> allPositions = new TreeSet<Integer>();

                for (CigarElement ce : read.getCigar().getCigarElements()) {
                    if (!ce.getOperator().equals(CigarOperator.H) && !ce.getOperator().equals(CigarOperator.S)) {
                        if (ce.getOperator().equals(CigarOperator.I)) {
                            insertionPositions.add(currentPosition);
                            allPositions.add(currentPosition);

                            numInsertions++;
                        } else if (ce.getOperator().equals(CigarOperator.D)) {
                            deletionPositions.add(currentPosition);
                            allPositions.add(currentPosition);

                            numDeletions++;
                        }

                        currentPosition += ce.getLength();
                    }
                }

                List<String> pieces = new ArrayList<String>();

                if (md != null) {
                    Matcher m = p.matcher(md);

                    currentPosition = 0;
                    while (m.find()) {
                        String piece = md.substring(m.start(), m.end());

                        if (!piece.contains("^")) {
                            if (piece.matches("[A-Z]")) {
                                numMismatches++;
                                currentPosition++;

                                mismatchPositions.add(currentPosition);
                                allPositions.add(currentPosition);
                            } else {
                                numMatches += Integer.parseInt(piece);
                                currentPosition += Integer.parseInt(piece);
                            }
                        }

                        pieces.add(piece);
                    }
                }

                incrementErrorsArray(mismatchPositions, mismatches, mismatchesTotal, read.getReadLength());
                incrementErrorsArray(insertionPositions, insertions, insertionsTotal, read.getReadLength());
                incrementErrorsArray(deletionPositions, deletions, deletionsTotal, read.getReadLength());
                incrementErrorsArray(allPositions, all, allTotal, read.getReadLength());

                Map<String, String> entry = new LinkedHashMap<String, String>();
                entry.put("numMatches", String.valueOf(numMatches));
                entry.put("numMismatches", String.valueOf(numMismatches));
                entry.put("numInsertions", String.valueOf(numInsertions));
                entry.put("numDeletions", String.valueOf(numDeletions));

                tw2.addEntry(entry);
            }
        }

        TableWriter tw = new TableWriter(out);
        for (int i = 0; i < lengthMax; i++) {
            Map<String, String> entry = new LinkedHashMap<String, String>();

            entry.put("index", String.valueOf(i));

            entry.put("mismatch", String.valueOf(mismatches[i]));
            entry.put("mismatchNorm", String.valueOf(mismatchesTotal[i]));
            entry.put("mismatchRate", String.valueOf((float) mismatches[i] / (float) mismatchesTotal[i]));

            entry.put("insertion", String.valueOf(insertions[i]));
            entry.put("insertionNorm", String.valueOf(insertionsTotal[i]));
            entry.put("insertionRate", String.valueOf((float) insertions[i] / (float) insertionsTotal[i]));

            entry.put("deletion", String.valueOf(deletions[i]));
            entry.put("deletionNorm", String.valueOf(deletionsTotal[i]));
            entry.put("deletionRate", String.valueOf((float) deletions[i] / (float) deletionsTotal[i]));

            entry.put("all", String.valueOf(all[i]));
            entry.put("allNorm", String.valueOf(allTotal[i]));
            entry.put("allRate", String.valueOf((float) all[i] / (float) allTotal[i]));

            if (allTotal[i] > 0) {
                tw.addEntry(entry);
            }
        }
    }
}
