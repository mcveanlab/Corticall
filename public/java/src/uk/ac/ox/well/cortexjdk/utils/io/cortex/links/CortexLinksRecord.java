package uk.ac.ox.well.cortexjdk.utils.io.cortex.links;

import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;

import java.util.ArrayList;
import java.util.List;

public class CortexLinksRecord {
    private String kmer;
    private List<CortexJunctionsRecord> cjs;

    public CortexLinksRecord(String skmer, List<CortexJunctionsRecord> cjs) {
        this.kmer = skmer;
        this.cjs = cjs;
    }

    public CortexLinksRecord(byte[] block) {
        String[] lines = new String(block).split("\\n");

        String[] kmerLine = lines[0].split("\\s+");

        kmer = kmerLine[0];
        int numLinks = Integer.valueOf(kmerLine[1]);

        cjs = new ArrayList<>();
        for (int i = 0; i < numLinks; i++) {
            String[] linkLine = lines[i + 1].split("\\s+");

            String orientation = linkLine[0];
            int numKmers = Integer.valueOf(linkLine[1]);

            String[] covStrs = linkLine[2].split(",");
            int[] coverages = new int[covStrs.length];

            for (int c = 0; c < coverages.length; c++) {
                coverages[c] = Integer.valueOf(covStrs[c]);
            }

            String junctions = linkLine[3];

            CortexJunctionsRecord cj = new CortexJunctionsRecord(orientation.equals("F"), numKmers, junctions.length(), coverages, junctions);
            cjs.add(cj);
        }
    }

    public String toString() {
        StringBuilder record = new StringBuilder();

        record.append(kmer).append(" ").append(cjs.size()).append("\n");

        int i = 0;
        for (CortexJunctionsRecord cj : cjs) {
            i++;

            record.append(cj);

            if (i < cjs.size()) { record.append("\n"); }
        }

        return record.toString();
    }

    public CortexByteKmer getKmerAsByteKmer() { return new CortexByteKmer(kmer.getBytes()); }
    public String getKmerAsString() { return kmer; }
    public CortexKmer getKmer() { return new CortexKmer(kmer); }
    public List<CortexJunctionsRecord> getJunctions() { return cjs; }

    public boolean equals(Object obj) {
        if (obj instanceof CortexLinksRecord) {
            List<CortexJunctionsRecord> js = getJunctions();

            boolean junctionsAreEqual = true;

            if (cjs.size() == js.size()) {
                for (int i = 0; i < js.size(); i++) {
                    if (!cjs.get(i).equals(js.get(i))) {
                        junctionsAreEqual = false;
                        break;
                    }
                }
            }

            return (junctionsAreEqual && kmer.equals(((CortexLinksRecord) obj).getKmerAsString()));
        }

        return false;
    }
}
