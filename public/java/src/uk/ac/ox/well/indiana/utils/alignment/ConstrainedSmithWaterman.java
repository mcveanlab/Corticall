package uk.ac.ox.well.indiana.utils.alignment;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.util.ArrayList;
import java.util.List;

public class ConstrainedSmithWaterman {
    private String reference;
    private String query;
    private int match;
    private int mismatch;
    private int gap;

    private int[][] matrix;

    private int alignmentScore = -1;
    private int alignmentStart;
    private String refAligned;
    private String queryAligned;

    public ConstrainedSmithWaterman(String reference, String query) {
        initialize(reference, query, 1, -1, -1);
    }

    public ConstrainedSmithWaterman(String reference, String query, int match, int mismatch, int gap) {
        initialize(reference, query, match, mismatch, gap);
    }

    private void initialize(String reference, String query, int match, int mismatch, int gap) {
        this.reference = reference;
        this.query = query;
        this.match = match;
        this.mismatch = mismatch;
        this.gap = gap;

        this.matrix = new int[query.length() + 1][reference.length() + 1];
    }

    public void addConstrainingKmer(String kmer) {
        int refStart = reference.indexOf(kmer) + 1;
        int queryStart = query.indexOf(kmer) + 1;

        if (refStart < 0 || queryStart < 0) {
            String rc = SequenceUtils.getReverseComplement(kmer);

            if (refStart < 0) {
                refStart = reference.indexOf(rc) + 1;
            }

            if (queryStart < 0) {
                queryStart = query.indexOf(rc) + 1;
            }
        }

        if (refStart >= 0 && queryStart >= 0) {
            for (int r = refStart, q = queryStart; r < kmer.length(); r++, q++) {
                matrix[q][r] = -1;
            }
        }
    }

    private int fillTable(int q, int r) {
        int value;

        if (q == 0 || r == 0) {
            value = 0;
        } else if (this.matrix[q][r] == -1) {
            value = fillTable(r - 1, q - 1) + match;
        } else if (this.matrix[q][r] > 0) {
            value = this.matrix[q][r];
        } else {
            int matchMismatchScore = fillTable(q - 1, r - 1) + (query.charAt(q - 1) == reference.charAt(r - 1) ? match : mismatch);
            int deletionScore = fillTable(q - 1, r) + gap;
            int insertionScore = fillTable(q, r - 1) + gap;

            value = Math.max(0, Math.max(matchMismatchScore, Math.max(insertionScore, deletionScore)));
        }

        //System.out.println("q=" + q + " r=" + r + " value=" + value);

        this.matrix[q][r] = value;

        return value;
    }

    private void walkTable() {
        fillTable(matrix.length - 1, matrix[0].length - 1);

        /*
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
        */

        StringBuilder queryAligned = new StringBuilder();
        StringBuilder refAligned = new StringBuilder();

        int maxScore = 0;
        int maxScoreR = 0;
        int maxScoreQ = 0;
        for (int q = 0; q < matrix.length; q++) {
            for (int r = 0; r < matrix[0].length; r++) {
                int cellScore = matrix[q][r];

                if (cellScore > maxScore) {
                    maxScore = cellScore;
                    maxScoreQ = q;
                    maxScoreR = r;
                }
            }
        }

        int q = maxScoreQ;
        int r = maxScoreR;

        queryAligned.insert(0, query.charAt(q - 1));
        refAligned.insert(0, reference.charAt(r - 1));

        while (q > 0 && r > 0) {
            int matchScore = matrix[q - 1][r - 1];
            int insertionScore = matrix[q][r - 1];
            int deletionScore = matrix[q - 1][r];

            if (insertionScore > matchScore && insertionScore > deletionScore) {
                r = r - 1;

                if (q > 0 && r > 0) {
                    queryAligned.insert(0, '-');
                    refAligned.insert(0, reference.charAt(r - 1));
                }
            } else if (deletionScore > insertionScore && deletionScore > matchScore) {
                q = q - 1;

                if (q > 0 && r > 0) {
                    queryAligned.insert(0, query.charAt(q - 1));
                    refAligned.insert(0, '-');
                }
            } else {
                q = q - 1;
                r = r - 1;

                if (q > 0 && r > 0) {
                    queryAligned.insert(0, query.charAt(q - 1));
                    refAligned.insert(0, reference.charAt(r - 1));
                }
            }
        }

        this.alignmentStart = r;
        this.queryAligned = queryAligned.toString();
        this.refAligned = refAligned.toString();
    }

    public int[][] getScoreTable() {
        if (this.alignmentScore == -1) {
            walkTable();
        }

        return this.matrix;
    }

    public int getAlignmentScore() {
        if (this.alignmentScore == -1) {
            walkTable();
        }

        return this.alignmentScore;
    }

    public int getAlignmentStart() {
        if (this.alignmentScore == -1) {
            walkTable();
        }

        return this.alignmentStart;
    }

    public String[] getAlignment() {
        if (this.alignmentScore == -1) {
            walkTable();
        }

        return new String[] { refAligned, queryAligned };
    }

    public Cigar getCigar() {
        List<CigarElement> cigarElements = new ArrayList<CigarElement>();

        String[] aligned = getAlignment();

        int initialSoftClipLength = query.indexOf(aligned[1].replaceAll("-", ""));
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

        int remainingSoftClipLength = query.length() - (initialSoftClipLength + aligned[1].replaceAll("-", "").length());
        if (remainingSoftClipLength > 0) {
            cigarElements.add(new CigarElement(remainingSoftClipLength, CigarOperator.SOFT_CLIP));
        }

        return new Cigar(cigarElements);
    }
}
