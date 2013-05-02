package uk.ac.ox.well.indiana.utils.assembly;

import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexMap;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.util.*;

public class CortexGraphWalker2 {
    private CortexMap cortexMap;

    public CortexGraphWalker2(CortexMap cortexMap) {
        this.cortexMap = cortexMap;
    }

    public CortexKmer buildContig(int color, CortexKmer panelKmer) {
        if (cortexMap.containsKey(panelKmer)) {
            StringBuilder contig = new StringBuilder();

            Set<CortexKmer> seenKmers = new HashSet<CortexKmer>();

            // start with left edges
            CortexRecord record = cortexMap.get(panelKmer);
            Collection<Byte> leftEdges = record.getInEdgesAsBytes(color);
            byte[] currentKmer = record.getKmerAsBytes();

            contig.append(new String(currentKmer));

            while (record != null && !seenKmers.contains(record.getKmer()) && leftEdges.size() == 1) {
                seenKmers.add(record.getKmer());

                byte[] newKmer = new byte[currentKmer.length];

                newKmer[0] = leftEdges.iterator().next();
                System.arraycopy(currentKmer, 0, newKmer, 1, currentKmer.length - 1);

                contig.insert(0, SequenceUtils.nucleotideByteToString(newKmer[0]));

                CortexKmer leftKmer = new CortexKmer(newKmer);
                record = cortexMap.get(leftKmer);

                if (record != null) {
                    leftEdges = leftKmer.isFlipped() ? record.getOutEdgesComplementAsBytes(color) : record.getInEdgesAsBytes(color);
                    currentKmer = newKmer;
                }
            }

            // now do right edges
            record = cortexMap.get(panelKmer);
            Collection<Byte> rightEdges = record.getOutEdgesAsBytes(color);
            currentKmer = record.getKmerAsBytes();

            seenKmers.remove(panelKmer);

            while (record != null && !seenKmers.contains(record.getKmer()) && rightEdges.size() == 1) {
                seenKmers.add(record.getKmer());

                byte[] newKmer = new byte[currentKmer.length];

                System.arraycopy(currentKmer, 1, newKmer, 0, currentKmer.length - 1);
                newKmer[currentKmer.length - 1] = rightEdges.iterator().next();

                contig.append(SequenceUtils.nucleotideByteToString(newKmer[currentKmer.length - 1]));

                CortexKmer rightKmer = new CortexKmer(newKmer);
                record = cortexMap.get(rightKmer);

                if (record != null) {
                    rightEdges = rightKmer.isFlipped() ? record.getInEdgesComplementAsBytes(color) : record.getOutEdgesAsBytes(color);
                    currentKmer = newKmer;
                }
            }

            return new CortexKmer(contig.toString());
        }

        return null;
    }

    public Map<CortexKmer, Set<CortexKmer>> buildContigs(int color, Set<CortexKmer> panelKmers) {
        Set<CortexKmer> usedKmers = new HashSet<CortexKmer>();

        Map<CortexKmer, Set<CortexKmer>> contigs = new HashMap<CortexKmer, Set<CortexKmer>>();

        for (CortexKmer panelKmer : panelKmers) {
            if (!usedKmers.contains(panelKmer)) {
                CortexKmer contig = buildContig(color, panelKmer);

                if (contig != null) {
                    Set<CortexKmer> seedKmers = new HashSet<CortexKmer>();

                    byte[] contigBytes = contig.getKmerAsBytes();

                    for (int i = 0; i <= contigBytes.length - panelKmer.length(); i++) {
                        byte[] kmer = new byte[panelKmer.length()];
                        System.arraycopy(contigBytes, i, kmer, 0, panelKmer.length());

                        CortexKmer ckmer = new CortexKmer(kmer);
                        if (panelKmers.contains(ckmer)) {
                            seedKmers.add(ckmer);
                        }
                    }

                    usedKmers.addAll(seedKmers);

                    contigs.put(contig, seedKmers);
                }
            }
        }

        return contigs;
    }
}
