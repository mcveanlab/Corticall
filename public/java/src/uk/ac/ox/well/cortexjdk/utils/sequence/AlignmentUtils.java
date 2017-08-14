package uk.ac.ox.well.cortexjdk.utils.sequence;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class AlignmentUtils {
    private AlignmentUtils() {}

    private static Pattern p = Pattern.compile("[0-9]+|[A-Z]|\\^[A-Z]+");

    private static int computeInitialPosition(SAMRecord read) {
        int currentPosition = 0;

        if (read.getCigar().getCigarElements().size() > 0 && read.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.S)) {
            currentPosition = read.getCigar().getCigarElement(0).getLength();
        }

        return currentPosition;
    }

    public static Set<Integer> getInsertionPositions(SAMRecord read) {
        int currentPosition = computeInitialPosition(read);

        Set<Integer> insertionPositions = new TreeSet<>();

        for (CigarElement ce : read.getCigar().getCigarElements()) {
            if (!ce.getOperator().equals(CigarOperator.H) && !ce.getOperator().equals(CigarOperator.S)) {
                if (ce.getOperator().equals(CigarOperator.I)) {
                    insertionPositions.add(currentPosition);
                }

                currentPosition += ce.getLength();
            }
        }

        return insertionPositions;
    }

    public static Set<Integer> getDeletionPositions(SAMRecord read) {
        int currentPosition = computeInitialPosition(read);

        Set<Integer> deletionPositions = new TreeSet<>();

        for (CigarElement ce : read.getCigar().getCigarElements()) {
            if (!ce.getOperator().equals(CigarOperator.H) && !ce.getOperator().equals(CigarOperator.S)) {
                if (ce.getOperator().equals(CigarOperator.D)) {
                    deletionPositions.add(currentPosition);
                }

                currentPosition += ce.getLength();
            }
        }

        return deletionPositions;
    }

    public static Set<Integer> getMismatchPositions(SAMRecord read) {
        Set<Integer> mismatchPositions = new TreeSet<>();

        String md = read.getStringAttribute("MD");
        if (md != null) {
            Matcher m = p.matcher(md);

            int currentPosition = computeInitialPosition(read);

            if (read.getCigar().getCigarElements().size() > 0 && read.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.S)) {
                currentPosition = read.getCigar().getCigarElement(0).getLength();
            }

            while (m.find()) {
                String piece = md.substring(m.start(), m.end());

                if (!piece.contains("^")) {
                    if (piece.matches("[A-Z]")) {
                        mismatchPositions.add(currentPosition);
                        currentPosition++;
                    } else {
                        currentPosition += Integer.parseInt(piece);
                    }
                }
            }
        }

        return mismatchPositions;
    }

    public static Set<Integer> getDifferencesPositions(SAMRecord read) {
        Set<Integer> diffPositions = new TreeSet<>();

        diffPositions.addAll(getMismatchPositions(read));
        diffPositions.addAll(getInsertionPositions(read));
        diffPositions.addAll(getDeletionPositions(read));

        return diffPositions;
    }

    public static List<Integer> getInsertionSizes(SAMRecord read) {
        List<Integer> insertionSizes = new ArrayList<>();

        for (CigarElement ce : read.getCigar().getCigarElements()) {
            if (ce.getOperator().equals(CigarOperator.I)) {
                insertionSizes.add(ce.getLength());
            }
        }

        return insertionSizes;
    }

    public static List<Integer> getDeletionSizes(SAMRecord read) {
        List<Integer> deletionSizes = new ArrayList<>();

        for (CigarElement ce : read.getCigar().getCigarElements()) {
            if (ce.getOperator().equals(CigarOperator.D)) {
                deletionSizes.add(ce.getLength());
            }
        }

        return deletionSizes;
    }

    public static List<CigarElement> getForwardCigar(SAMRecord alignment) {
        List<CigarElement> ces = new ArrayList<>();

        if (alignment.getReadNegativeStrandFlag()) {
            for (int i = alignment.getCigar().numCigarElements() - 1; i >= 0; i--) {
                ces.add(alignment.getCigar().getCigarElement(i));
            }
        } else {
            ces = alignment.getCigar().getCigarElements();
        }

        return ces;
    }
}
