package uk.ac.ox.well.cortexjdk.utils.io.cortex;

import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.*;

import java.io.File;
import java.util.Iterator;
import java.util.List;

/**
 * Created by kiran on 20/08/2017.
 */
public interface DeBruijnGraph extends Iterable<CortexRecord>, Iterator<CortexRecord> {
    // Seeking
    long position();
    void position(long i);

    // Iterating
    Iterator<CortexRecord> iterator();
    boolean hasNext();
    CortexRecord next();
    void remove();
    void close();

    // Finding specific records
    CortexRecord getRecord(long i);
    CortexRecord findRecord(byte[] bk);
    CortexRecord findRecord(CortexByteKmer bk);
    CortexRecord findRecord(CortexKmer ck);
    CortexRecord findRecord(String sk);

    // Graph information
    File getFile();
    CortexHeader getHeader();
    int getVersion();
    int getKmerSize();
    int getKmerBits();
    int getNumColors();
    long getNumRecords();

    // Color information
    List<CortexColor> getColors();
    boolean hasColor(int color);
    CortexColor getColor(int color);
    int getColorForSampleName(String sampleName);
    List<Integer> getColorsForSampleNames(List<String> sampleNames);
    String getSampleName(int color);

    String toString();
}

