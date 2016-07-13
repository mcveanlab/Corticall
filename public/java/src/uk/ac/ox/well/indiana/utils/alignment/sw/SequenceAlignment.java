package uk.ac.ox.well.indiana.utils.alignment.sw;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Paul Reiners
 *
 */
public abstract class SequenceAlignment extends DynamicProgramming {
    protected int match;
    protected int mismatch;
    protected int space;
    protected String[] alignments;

    public SequenceAlignment(String sequence1, String sequence2) {
        this(sequence1, sequence2, 1, -1, -1);
    }

    public SequenceAlignment(String sequence1, String sequence2, int match, int mismatch, int gap) {
        super(sequence1, sequence2);

        this.match = match;
        this.mismatch = mismatch;
        this.space = gap;
    }

    protected Object getTraceback() {
        StringBuilder align1Buf = new StringBuilder();
        StringBuilder align2Buf = new StringBuilder();
        Cell currentCell = getTracebackStartingCell();

        while (traceBackIsNotDone(currentCell)) {
            if (currentCell.getRow() - currentCell.getPrevCell().getRow() == 1) {
                align2Buf.insert(0, sequence2.charAt(currentCell.getRow() - 1));
            } else {
                align2Buf.insert(0, '-');
            }
            if (currentCell.getCol() - currentCell.getPrevCell().getCol() == 1) {
                align1Buf.insert(0, sequence1.charAt(currentCell.getCol() - 1));
            } else {
                align1Buf.insert(0, '-');
            }
            currentCell = currentCell.getPrevCell();
        }

        String[] alignments = new String[] { align1Buf.toString(), align2Buf.toString() };

        return alignments;
    }

    protected abstract boolean traceBackIsNotDone(Cell currentCell);

    public int getAlignmentScore() {
        if (alignments == null) {
            getAlignment();
        }

        int score = 0;
        for (int i = 0; i < alignments[0].length(); i++) {
            char c1 = alignments[0].charAt(i);
            char c2 = alignments[1].charAt(i);
            if (c1 == '-' || c2 == '-') {
                score += space;
            } else if (c1 == c2) {
                score += match;
            } else {
                score += mismatch;
            }
        }

        return score;
    }

    public String[] getAlignment() {
        ensureTableIsFilledIn();
        alignments = (String[]) getTraceback();
        return alignments;
    }

    public int getAlignmentStart() {
        String[] aligned = getAlignment();

        int initialSoftClipLength = sequence2.indexOf(aligned[1].replaceAll("-", ""));

        int alignmentStart = sequence1.indexOf(aligned[0].replaceAll("-", "")) - sequence2.indexOf(aligned[1].replaceAll("-", "")) + initialSoftClipLength;

        return alignmentStart;
    }

    public Cigar getCigar() {
        List<CigarElement> cigarElements = new ArrayList<>();

        String[] aligned = getAlignment();

        int initialSoftClipLength = sequence2.indexOf(aligned[1].replaceAll("-", ""));
        if (initialSoftClipLength > 0) {
            cigarElements.add(new CigarElement(initialSoftClipLength, CigarOperator.SOFT_CLIP));
        }

        int currentElementLength = 0;
        CigarOperator currentElementOperator = null;

        for (int i = 0; i < aligned[0].length(); i++) {
            CigarOperator nextElementOperator;

            if (aligned[0].charAt(i) == '-') {
                nextElementOperator = CigarOperator.INSERTION;
            } else if (aligned[1].charAt(i) == '-') {
                nextElementOperator = CigarOperator.DELETION;
            } else {
                nextElementOperator = CigarOperator.MATCH_OR_MISMATCH;
            }

            if (currentElementOperator == null) {
                currentElementOperator = nextElementOperator;
                currentElementLength = 0;
            } else if (currentElementOperator != nextElementOperator) {
                cigarElements.add(new CigarElement(currentElementLength, currentElementOperator));

                currentElementOperator = nextElementOperator;
                currentElementLength = 0;
            }

            currentElementLength++;
        }

        cigarElements.add(new CigarElement(currentElementLength, currentElementOperator));

        int remainingSoftClipLength = sequence2.length() - (initialSoftClipLength + aligned[1].replaceAll("-", "").length());
        if (remainingSoftClipLength > 0) {
            cigarElements.add(new CigarElement(remainingSoftClipLength, CigarOperator.SOFT_CLIP));
        }

        return new Cigar(cigarElements);
    }

    protected abstract Cell getTracebackStartingCell();
}

