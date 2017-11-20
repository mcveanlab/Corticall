package uk.ac.ox.well.cortexjdk.utils.io.graph.links;

import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexColor;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexHeader;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexBinaryKmer;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

public class CortexLinksMap implements uk.ac.ox.well.cortexjdk.utils.io.graph.ConnectivityAnnotations {
    private CortexLinksIterable cortexGraphLinks;
    private Map<CortexBinaryKmer, CortexLinksRecord> recordHash;
    private CortexHeader header;

    public CortexLinksMap(String cortexLinksPath) { initialize(new File(cortexLinksPath)); }

    public CortexLinksMap(File cortexLinksFile) { initialize(cortexLinksFile); }

    private void initialize(File cortexLinksFile) {
        this.cortexGraphLinks = new CortexLinksIterable(cortexLinksFile);

        header = new CortexHeader();
        header.setKmerSize(cortexGraphLinks.getKmerSize());
        header.setKmerBits(CortexRecord.getKmerBits(cortexGraphLinks.getKmerSize()));
        header.setVersion(cortexGraphLinks.getVersion());
        header.setNumColors(cortexGraphLinks.getNumColors());

        for (int c = 0; c < header.getNumColors(); c++) {
            CortexColor cc = new CortexColor();
            cc.setSampleName(cortexGraphLinks.getColor(c).getSampleName());
        }

        recordHash = new HashMap<>((int) cortexGraphLinks.getNumLinks());

        for (CortexLinksRecord clr : cortexGraphLinks) {
            recordHash.put(new CortexBinaryKmer(clr.getKmer().getKmerAsBytes()), clr);
        }
    }

    @Override
    public int size() { return recordHash.size(); }

    @Override
    public boolean isEmpty() { return recordHash.isEmpty(); }

    @Override
    public boolean containsKey(Object key) { return recordHash.containsKey(convert(key)); }

    @Override
    public CortexLinksRecord get(Object key) { return recordHash.get(convert(key)); }

    @Override
    public CortexHeader getHeader() { return header; }
}
