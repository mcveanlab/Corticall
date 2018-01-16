package uk.ac.ox.well.cortexjdk.utils.io.graph.cortex;

import uk.ac.ox.well.cortexjdk.Main;
import uk.ac.ox.well.cortexjdk.utils.io.graph.DeBruijnGraph;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexBinaryKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;

import java.io.File;
import java.util.*;

public class CortexMap implements DeBruijnGraph {
    private CortexGraph graph;

    private Map<CortexBinaryKmer, CortexRecord> recs = new HashMap<>();

    public CortexMap(String cortexFilePath) { loadGraph(new File(cortexFilePath)); }
    public CortexMap(File cortexFile) { loadGraph(cortexFile); }

    private void loadGraph(File cortexFile) {
        this.graph = new CortexGraph(cortexFile);

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Pre-loading graph...")
                .message("records loaded")
                .updateRecord(10000000)
                .make(Main.getLogger());

        for (CortexRecord cr : graph) {
            recs.put(cr.getCortexBinaryKmer(), cr);

            pm.update();
        }
    }

    @Override
    public long position() {
        return graph.position();
    }

    @Override
    public void position(long i) {
        graph.position(i);
    }

    @Override
    public Iterator<CortexRecord> iterator() {
        return graph.iterator();
    }

    @Override
    public boolean hasNext() {
        return graph.hasNext();
    }

    @Override
    public CortexRecord next() {
        return graph.next();
    }

    @Override
    public void remove() {
        graph.remove();
    }

    @Override
    public void close() {
        graph.close();
    }

    @Override
    public CortexRecord getRecord(long i) {
        return graph.getRecord(i);
    }

    private CortexRecord get(CortexBinaryKmer cbk) {
        return recs.getOrDefault(cbk, null);
    }

    @Override
    public CortexRecord findRecord(byte[] bk) {
        return get(new CortexBinaryKmer(bk));
    }

    @Override
    public CortexRecord findRecord(CortexByteKmer bk) {
        return get(new CortexBinaryKmer(bk.getKmer()));
    }

    @Override
    public CortexRecord findRecord(CanonicalKmer ck) {
        return get(new CortexBinaryKmer(ck.getKmerAsBytes()));
    }

    @Override
    public CortexRecord findRecord(String sk) {
        return get(new CortexBinaryKmer(sk.getBytes()));
    }

    @Override
    public File getFile() {
        return graph.getFile();
    }

    @Override
    public CortexHeader getHeader() {
        return graph.getHeader();
    }

    @Override
    public int getVersion() {
        return graph.getVersion();
    }

    @Override
    public int getKmerSize() {
        return graph.getKmerSize();
    }

    @Override
    public int getKmerBits() {
        return graph.getKmerBits();
    }

    @Override
    public int getNumColors() {
        return graph.getNumColors();
    }

    @Override
    public long getNumRecords() {
        return graph.getNumRecords();
    }

    @Override
    public List<CortexColor> getColors() {
        return graph.getColors();
    }

    @Override
    public boolean hasColor(int color) {
        return graph.hasColor(color);
    }

    @Override
    public CortexColor getColor(int color) {
        return graph.getColor(color);
    }

    @Override
    public int getColorForSampleName(String sampleName) {
        return graph.getColorForSampleName(sampleName);
    }

    @Override
    public List<Integer> getColorsForSampleNames(Collection<String> sampleNames) {
        return graph.getColorsForSampleNames(sampleNames);
    }

    @Override
    public String getSampleName(int color) {
        return graph.getSampleName(color);
    }
}
