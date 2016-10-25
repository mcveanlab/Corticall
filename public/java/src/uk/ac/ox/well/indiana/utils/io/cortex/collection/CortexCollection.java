package uk.ac.ox.well.indiana.utils.io.cortex.collection;

import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

public class CortexCollection implements Iterable<CortexRecord>, Iterator<CortexRecord> {
    private CortexRecord[] nextRecs;

    public CortexCollection(String collectionFileString) {
    }

    public CortexGraph getGraph(int color) {
        return null;
    }

    public int getNumColors() {
        return 0;
    }
    public int getKmerSize() {
        return 0;
    }
    public int getKmerBits() {
        return 0;
    }

    public CortexRecord findRecord(String kmer) {
        return null;
    }

    private void moveToBeginningOfRecordsSection() {
    }

    @Override
    public Iterator<CortexRecord> iterator() {
        moveToBeginningOfRecordsSection();

        return this;
    }

    @Override
    public boolean hasNext() {
        return false;
    }

    @Override
    public CortexRecord next() {
        return null;
    }
}
