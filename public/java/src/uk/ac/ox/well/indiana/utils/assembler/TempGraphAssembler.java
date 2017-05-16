package uk.ac.ox.well.indiana.utils.assembler;

import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.*;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Created by kiran on 13/05/2017.
 */
public class TempGraphAssembler {
    private TempGraphAssembler() {}

    public static CortexGraph buildGraph(Map<String, Collection<String>> haplotypeLists, int kmerSize) {
        File tempFile = null;
        try {
            tempFile = File.createTempFile("tempgraph", ".ctx");
        } catch (IOException e) {
            throw new IndianaException("Could not get a temp file for graph creation");
        }
        tempFile.deleteOnExit();

        CortexGraphWriter cgw = new CortexGraphWriter(tempFile);
        cgw.setHeader(constructCortexHeader(haplotypeLists.keySet(), kmerSize));

        Map<CortexKmer, CortexRecord> crs = new TreeMap<>();

        int c = 0;
        for (String sampleName : haplotypeLists.keySet()) {
            for (String sequence : haplotypeLists.get(sampleName)) {
                sequence = sequence.toUpperCase();

                for (int i = 0; i <= sequence.length() - kmerSize; i++) {
                    String sk = sequence.substring(i, i + kmerSize);

                    String prevBase = i == 0 ? null : sequence.substring(i - 1, i);
                    String nextBase = i == sequence.length() - kmerSize ? null : sequence.substring(i + kmerSize, i + kmerSize + 1);

                    updateRecord(crs, haplotypeLists.size(), c, sk, prevBase, nextBase);
                }
            }

            c++;
        }

        for (CortexKmer cr : crs.keySet()) {
            cgw.addRecord(crs.get(cr));
        }

        cgw.close();

        return new CortexGraph(tempFile);
    }

    public static void updateRecord(Map<CortexKmer, CortexRecord> crs, int numColors, int color, String sk, String prevBase, String nextBase) {
        CortexKmer ck = new CortexKmer(sk);

        List<Integer> coverageList = new ArrayList<>();
        List<Set<String>> inEdgesList = new ArrayList<>();
        List<Set<String>> outEdgesList = new ArrayList<>();

        CortexRecord oldRc = crs.containsKey(ck) ? crs.get(ck) : null;

        for (int c = 0; c < numColors; c++) {
            int cov = 0;
            Set<String> inEdges = new HashSet<>();
            Set<String> outEdges = new HashSet<>();

            if (oldRc != null) {
                cov += oldRc.getCoverage(c);

                inEdges.addAll(oldRc.getInEdgesAsStrings(c, false));
                outEdges.addAll(oldRc.getOutEdgesAsStrings(c, false));
            }

            if (c == color) {
                cov++;

                if (!ck.isFlipped()) {
                    if (prevBase != null) { inEdges.add(prevBase); }
                    if (nextBase != null) { outEdges.add(nextBase); }
                } else {
                    if (nextBase != null) { inEdges.add(SequenceUtils.reverseComplement(nextBase)); }
                    if (prevBase != null) { outEdges.add(SequenceUtils.reverseComplement(prevBase)); }
                }
            }

            coverageList.add(cov);
            inEdgesList.add(inEdges);
            outEdgesList.add(outEdges);
        }

        crs.put(ck, new CortexRecord(SequenceUtils.alphanumericallyLowestOrientation(sk), coverageList, inEdgesList, outEdgesList));
    }

    @NotNull
    public static CortexHeader constructCortexHeader(Set<String> colors, int kmerSize) {
        CortexHeader ch = new CortexHeader();
        ch.setVersion(6);
        ch.setNumColors(colors.size());
        ch.setKmerSize(kmerSize);
        ch.setKmerBits(CortexRecord.getKmerBits(kmerSize));

        for (String sampleName : colors) {
            CortexColor col = new CortexColor();

            col.setSampleName(sampleName);
            col.setCleanedAgainstGraph(false);
            col.setCleanedAgainstGraphName("");
            col.setErrorRate(0.0);
            col.setLowCovgKmersRemoved(false);
            col.setLowCovgSupernodesRemoved(false);
            col.setLowCovKmerThreshold(0);
            col.setTipClippingApplied(false);
            col.setMeanReadLength(0);
            col.setTotalSequence(0);

            ch.addColor(col);
        }

        return ch;
    }

}
