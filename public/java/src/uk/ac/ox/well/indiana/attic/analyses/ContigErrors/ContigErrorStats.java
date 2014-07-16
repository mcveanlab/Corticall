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

public class ContigErrorStats extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (in BAM format)")
    public SAMFileReader CONTIGS;

    @Argument(fullName="lengthMin", shortName="lmin", doc="Minimum contig length to process")
    public Integer LENGTH_MIN;

    @Argument(fullName="lengthMax", shortName="lmax", doc="Maximum contig length to process")
    public Integer LENGTH_MAX;

    @Argument(fullName="countInBothOrientations", shortName="cb", doc="Count mismatches in contigs in both orientations of the contig")
    public Boolean COUNT_IN_BOTH_ORIENTATIONS = false;

    @Output
    public PrintStream out;

    @Output(fullName="out2", shortName="o2")
    public PrintStream out2;

    private Set<Integer> reverseErrorSet(Set<Integer> errorSet, int contigLength) {
        Set<Integer> newErrorSet = new TreeSet<Integer>();
        for (Integer errorPos : errorSet) {
            newErrorSet.add(contigLength - errorPos);
        }

        return newErrorSet;
    }

    private void incrementErrorsArray(Set<Integer> errorSet, int[] errors, int[] norms, int contigLength) {
        for (Integer error : errorSet) {
            errors[error]++;
        }

        for (int i = 0; i < contigLength; i++) {
            norms[i]++;
        }
    }

    @Override
    public void execute() {
        Pattern p = Pattern.compile("[0-9]+|[A-Z]|\\^[A-Z]+");

        int[] mismatches = new int[LENGTH_MAX];
        int[] mismatchesTotal = new int[LENGTH_MAX];
        int[] insertions = new int[LENGTH_MAX];
        int[] insertionsTotal = new int[LENGTH_MAX];
        int[] deletions = new int[LENGTH_MAX];
        int[] deletionsTotal = new int[LENGTH_MAX];
        int[] all = new int[LENGTH_MAX];
        int[] allTotal = new int[LENGTH_MAX];

        TableWriter tw2 = new TableWriter(out2);

        for (SAMRecord read : CONTIGS) {
            if (read.getAlignmentStart() > 0 && read.getReadLength() >= LENGTH_MIN && read.getReadLength() < LENGTH_MAX) {
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

                //String joined = Joiner.on(",").join(pieces);
                //out.println(Joiner.on("\t").useForNull("").join(read.getCigarString(), md, joined, numMatches, numMismatches, numDeletions, numInsertions));
                //out.println("I " + Joiner.on(", ").join(insertionPositions));
                //out.println("D " + Joiner.on(", ").join(deletionPositions));
                //out.println("M " + Joiner.on(", ").join(mismatchPositions));

                incrementErrorsArray(mismatchPositions, mismatches, mismatchesTotal, read.getReadLength());
                incrementErrorsArray(insertionPositions, insertions, insertionsTotal, read.getReadLength());
                incrementErrorsArray(deletionPositions, deletions, deletionsTotal, read.getReadLength());
                incrementErrorsArray(allPositions, all, allTotal, read.getReadLength());

                if (COUNT_IN_BOTH_ORIENTATIONS) {
                    incrementErrorsArray(reverseErrorSet(mismatchPositions, read.getReadLength()), mismatches, mismatchesTotal, read.getReadLength());
                    incrementErrorsArray(reverseErrorSet(insertionPositions, read.getReadLength()), insertions, insertionsTotal, read.getReadLength());
                    incrementErrorsArray(reverseErrorSet(deletionPositions, read.getReadLength()), deletions, deletionsTotal, read.getReadLength());
                    incrementErrorsArray(reverseErrorSet(allPositions, read.getReadLength()), all, allTotal, read.getReadLength());
                }

                Map<String, String> entry = new LinkedHashMap<String, String>();
                entry.put("numMatches", String.valueOf(numMatches));
                entry.put("numMismatches", String.valueOf(numMismatches));
                entry.put("numInsertions", String.valueOf(numInsertions));
                entry.put("numDeletions", String.valueOf(numDeletions));

                tw2.addEntry(entry);
            }
        }

        TableWriter tw = new TableWriter(out);
        for (int i = 0; i < LENGTH_MAX; i++) {
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

            tw.addEntry(entry);
        }
    }
}
