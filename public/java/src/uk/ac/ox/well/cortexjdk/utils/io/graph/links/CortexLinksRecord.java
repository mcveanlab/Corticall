package uk.ac.ox.well.cortexjdk.utils.io.graph.links;

import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;

import java.util.*;

public class CortexLinksRecord {
    private String kmer;
    private Set<CortexJunctionsRecord> cjs = new HashSet<>();

    public CortexLinksRecord(String skmer, Collection<CortexJunctionsRecord> cjs) {
        this.kmer = skmer;
        this.cjs.addAll(cjs);
    }

    public CortexLinksRecord(byte[] block) {
        String[] lines = new String(block).split("\\n");

        String[] kmerLine = lines[0].split("\\s+");

        kmer = kmerLine[0];
        int numLinks = Integer.valueOf(kmerLine[1]);

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
    public CanonicalKmer getKmer() { return new CanonicalKmer(kmer); }
    public Collection<CortexJunctionsRecord> getJunctions() { return cjs; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        CortexLinksRecord that = (CortexLinksRecord) o;

        if (kmer != null ? !kmer.equals(that.kmer) : that.kmer != null) return false;
        return cjs != null ? cjs.equals(that.cjs) : that.cjs == null;
    }

    @Override
    public int hashCode() {
        int result = kmer != null ? kmer.hashCode() : 0;
        result = 31 * result + (cjs != null ? cjs.hashCode() : 0);
        return result;
    }
}
