package uk.ac.ox.well.indiana.utils.io.cortex;

import java.io.File;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import ch.qos.logback.classic.Logger;

public class CortexMap implements Map<CortexKmer, CortexRecord> {
    private CortexGraph cortexGraph;
    private Logger log;

    private Map<CortexKmer, CortexRecord> recordHash;

    public CortexMap(String cortexGraphPath) {
        initialize(new CortexGraph(cortexGraphPath), null);
    }

    public CortexMap(String cortexGraphPath, Logger log) {
        initialize(new CortexGraph(cortexGraphPath), log);
    }

    public CortexMap(File cortexGraphFile) {
        initialize(new CortexGraph(cortexGraphFile), null);
    }

    public CortexMap(File cortexGraphFile, Logger log) {
        initialize(new CortexGraph(cortexGraphFile), log);
    }

    public CortexMap(CortexGraph cortexGraph) {
        initialize(this.cortexGraph = cortexGraph, null);
    }

    public CortexMap(CortexGraph cortexGraph, Logger log) {
        initialize(this.cortexGraph = cortexGraph, log);
    }

    private void initialize(CortexGraph cortexGraph, Logger log) {
        this.cortexGraph = cortexGraph;
        this.log = log;

        recordHash = new HashMap<CortexKmer, CortexRecord>((int) this.cortexGraph.getNumRecords());

        if (log != null) { log.info("Loading Cortex records"); }

        int i = 0;
        for (CortexRecord cr : cortexGraph) {
            put(new CortexKmer(cr.getKmerAsBytes(), true), cr);

            if (log != null && i % (cortexGraph.getNumRecords() / 5) == 0) {
                log.info("Loaded {}/{} Cortex records", i + 1, cortexGraph.getNumRecords());
            }
            i++;
        }

        if (log != null) { log.info("Finished loading Cortex records"); }
    }

    public boolean containsKey(String kmer) {
        return containsKey(new CortexKmer(kmer));
    }

    public boolean containsKey(byte[] kmer) {
        return containsKey(new CortexKmer(kmer));
    }

    public CortexRecord get(String kmer) {
        return get(new CortexKmer(kmer));
    }

    public CortexRecord get(byte[] kmer) {
        return get(new CortexKmer(kmer));
    }

    @Override
    public int size() {
        return recordHash.size();
    }

    @Override
    public boolean isEmpty() {
        return recordHash.isEmpty();
    }

    @Override
    public boolean containsKey(Object key) {
        return recordHash.containsKey(key);
    }

    @Override
    public boolean containsValue(Object value) {
        return recordHash.containsValue(value);
    }

    @Override
    public CortexRecord get(Object key) {
        return recordHash.get(key);
    }

    @Override
    public CortexRecord put(CortexKmer key, CortexRecord value) {
        return recordHash.put(key, value);
    }

    @Override
    public CortexRecord remove(Object key) {
        return recordHash.remove(key);
    }

    @Override
    public void putAll(Map<? extends CortexKmer, ? extends CortexRecord> m) {
        recordHash.putAll(m);
    }

    @Override
    public void clear() {
        recordHash.clear();
    }

    @Override
    public Set<CortexKmer> keySet() {
        return recordHash.keySet();
    }

    @Override
    public Collection<CortexRecord> values() {
        return recordHash.values();
    }

    @Override
    public Set<Entry<CortexKmer, CortexRecord>> entrySet() {
        return recordHash.entrySet();
    }
}
