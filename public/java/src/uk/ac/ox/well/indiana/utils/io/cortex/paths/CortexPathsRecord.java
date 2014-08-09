package uk.ac.ox.well.indiana.utils.io.cortex.paths;

import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;

import java.util.List;

public class CortexPathsRecord {
    private String kmer;
    private List<CortexJunctionsRecord> cjs;

    public CortexPathsRecord(String skmer, List<CortexJunctionsRecord> cjs) {
        this.kmer = skmer;
        this.cjs = cjs;
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

    public String getKmerAsString() { return kmer; }
    public CortexKmer getKmer() { return new CortexKmer(kmer); }
    public List<CortexJunctionsRecord> getJunctions() { return cjs; }

    public boolean equals(Object obj) {
        if (obj instanceof CortexPathsRecord) {
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

            return (junctionsAreEqual && kmer.equals(((CortexPathsRecord) obj).getKmerAsString()));
        }

        return false;
    }
}
