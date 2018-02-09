package uk.ac.ox.well.cortexjdk.utils.alignment.graph;

import htsjdk.samtools.fastq.FastqRecord;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;

import java.util.*;

public class SingleEndAlignmentInfo {
    private FastqRecord fq;
    private String[] sks;
    private CortexRecord[] crs;
    private int kmerSize;

    public SingleEndAlignmentInfo(FastqRecord fq, String[] sks, CortexRecord[] crs) {
        this.fq = fq;
        this.sks = sks;
        this.crs = crs;

        for (CortexRecord cr : crs) {
            if (cr != null) {
                kmerSize = cr.getKmerSize();
                break;
            }
        }
    }

    @Override
    public String toString() {
        return "SingleEndAlignmentInfo{" +
                "fq=" + fq +
                ", crs=" + Arrays.toString(crs) +
                '}';
    }

    public CortexRecord[] getRecords() { return crs; }

    public FastqRecord getFastqRecord() {
        return fq;
    }

    public boolean isSplitRead(int baselineColor, int comparisonColor) {
        for (int i = 0; i < crs.length; i++) {
            if (crs[i] != null && crs[i].getCoverage(baselineColor) > 0 && crs[i].getCoverage(comparisonColor) == 0) {
                int numMissingKmers = 1;

                int edgeConflictLeftIndex = i;
                boolean edgeConflictLeft = false;
                for (int j = i - 1; j >= 0; j--) {
                    if (crs[j] != null) {
                        if (crs[j].getCoverage(baselineColor) > 0 && crs[j].getCoverage(comparisonColor) == 0) {
                            numMissingKmers++;
                        } else if (crs[j].getEdges()[baselineColor] != crs[j].getEdges()[comparisonColor]) {
                            edgeConflictLeft = true;
                            edgeConflictLeftIndex = j;
                            break;
                        }
                    }
                }

                int edgeConflictRightIndex = i;
                boolean edgeConflictRight = false;
                for (int j = i + 1; j < crs.length; j++) {
                    if (crs[j] != null) {
                        if (crs[j].getCoverage(baselineColor) > 0 && crs[j].getCoverage(comparisonColor) == 0) {
                            numMissingKmers++;
                        } else if (crs[j].getEdges()[baselineColor] != crs[j].getEdges()[comparisonColor]) {
                            edgeConflictRight = true;
                            edgeConflictRightIndex = j;
                            break;
                        }
                    }
                }

                float missingPct = 100.0f * ((float) numMissingKmers / ((float) (edgeConflictRightIndex - edgeConflictLeftIndex - 1)));

                if (missingPct >= 80 && edgeConflictLeft && edgeConflictRight) {
                    return true;
                }
            }
        }

        return false;
    }

    public Pair<String, String> splitRead(int baselineColor, int comparisonColor) {
        for (int i = 0; i < crs.length; i++) {
            if (crs[i] != null && crs[i].getCoverage(baselineColor) > 0 && crs[i].getCoverage(comparisonColor) == 0) {
                int numMissingKmers = 1;

                int edgeConflictLeftIndex = i;
                boolean edgeConflictLeft = false;
                for (int j = i - 1; j >= 0; j--) {
                    if (crs[j] != null) {
                        if (crs[j].getCoverage(baselineColor) > 0 && crs[j].getCoverage(comparisonColor) == 0) {
                            numMissingKmers++;
                        } else if (crs[j].getEdges()[baselineColor] != crs[j].getEdges()[comparisonColor]) {
                            edgeConflictLeft = true;
                            edgeConflictLeftIndex = j;
                            break;
                        }
                    }
                }

                int edgeConflictRightIndex = i;
                boolean edgeConflictRight = false;
                for (int j = i + 1; j < crs.length; j++) {
                    if (crs[j] != null) {
                        if (crs[j].getCoverage(baselineColor) > 0 && crs[j].getCoverage(comparisonColor) == 0) {
                            numMissingKmers++;
                        } else if (crs[j].getEdges()[baselineColor] != crs[j].getEdges()[comparisonColor]) {
                            edgeConflictRight = true;
                            edgeConflictRightIndex = j;
                            break;
                        }
                    }
                }

                float missingPct = 100.0f * ((float) numMissingKmers / ((float) (edgeConflictRightIndex - edgeConflictLeftIndex - 1)));

                if (missingPct >= 80 && edgeConflictLeft && edgeConflictRight) {
                    StringBuilder end1 = new StringBuilder();

                    for (i = 0; i < edgeConflictLeftIndex; i++) {
                        String sk = fq.getReadString().substring(i, i + kmerSize);

                        if (end1.length() == 0) {
                            end1.append(sk);
                        } else {
                            end1.append(sk.substring(sk.length() - 1, sk.length()));
                        }
                    }

                    StringBuilder end2 = new StringBuilder();
                    for (i = edgeConflictRightIndex; i < crs.length; i++) {
                        String sk = fq.getReadString().substring(i, i + kmerSize);

                        if (end2.length() == 0) {
                            end2.append(sk);
                        } else {
                            end2.append(sk.substring(sk.length() - 1, sk.length()));
                        }
                    }

                    return new Pair<>(end1.toString(), end2.toString());
                }
            }
        }

        return null;
    }

    public List<String> getAvailableKmers(int color) {
        List<String> kmers = new ArrayList<>();

        for (int i = 0; i < sks.length; i++) {
            if (crs[i] != null && crs[i].getCoverage(color) > 0 && crs[i].getInDegree(color) == 1 && crs[i].getOutDegree(color) == 1) {
                kmers.add(sks[i]);
            }
        }

        return kmers;
    }
}
