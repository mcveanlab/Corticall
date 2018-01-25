package uk.ac.ox.well.cortexjdk.utils.io.graph.links;

import htsjdk.samtools.util.BlockCompressedInputStream;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexColor;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexHeader;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.utils.BinaryFile;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexBinaryKmer;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by kiran on 14/09/2017.
 */
public class CortexLinksRandomAccess implements uk.ac.ox.well.cortexjdk.utils.io.graph.ConnectivityAnnotations {
    private BlockCompressedInputStream bi;
    private Map<CortexBinaryKmer, Pair<Long, Integer>> index;
    private CortexHeader header;
    private String source;

    public CortexLinksRandomAccess(String cortexLinksPath) { initialize(new File(cortexLinksPath)); }

    public CortexLinksRandomAccess(File cortexLinksFile) { initialize(cortexLinksFile); }

    private void initialize(File cortexLinksFile) {
        File cortexLinksIndex = new File(cortexLinksFile.getAbsolutePath() + ".idx");

        try {
            bi = new BlockCompressedInputStream(cortexLinksFile.getAbsoluteFile());
            index = new HashMap<>();

            BinaryFile bf = new BinaryFile(cortexLinksIndex, "r");

            byte[] magicWordStart = new byte[6];
            bf.read(magicWordStart);

            header = new CortexHeader();
            header.setNumColors(bf.readInt());
            header.setKmerSize(bf.readInt());
            header.setKmerBits(CortexRecord.getKmerBits(header.getKmerSize()));
            long numKmersInGraph = bf.readLong();
            long numKmersWithLinks = bf.readLong();
            long numLinkBytes = bf.readLong();

            byte[] source = new byte[bf.readInt()];
            bf.read(source);
            this.source = new String(source);

            for (int c = 0; c < header.getNumColors(); c++) {
                CortexColor cc = new CortexColor();

                byte[] sn = new byte[bf.readInt()];
                bf.read(sn);

                cc.setSampleName(new String(sn));
                header.addColor(cc);
            }

            byte[] magicWordEnd = new byte[6];
            bf.read(magicWordEnd);

            if (!Arrays.equals(magicWordStart, magicWordEnd)) {
                throw new CortexJDKException("Error in decoding Cortex links index");
            }

            for (int i = 0; i < numKmersWithLinks; i++) {
                long[] kmer = new long[header.getKmerBits()];
                for (int j = 0; j < header.getKmerBits(); j++) {
                    kmer[j] = bf.readLong();
                }

                long pos = bf.readLong();
                int length = bf.readInt();

                index.put(new CortexBinaryKmer(kmer), new Pair<>(pos, length));
            }
        } catch (IOException e) {
            throw new CortexJDKException("IOException", e);
        }
    }

    @Override
    public int size() { return index.size(); }

    @Override
    public boolean isEmpty() { return index.isEmpty(); }

    @Override
    public boolean containsKey(Object key) { return index.containsKey(convert(key)); }

    @Override
    public CortexLinksRecord get(Object key) {
        try {
            Pair<Long, Integer> p = index.get(convert(key));
            byte[] recbuf = new byte[p.getSecond()];

            bi.seek(p.getFirst());
            bi.read(recbuf);

            return new CortexLinksRecord(recbuf);
        } catch (IOException e) {
            throw new CortexJDKException("Failed to load links record from disk", e);
        }
    }

    @Override
    public CortexHeader getHeader() { return header; }

    @Override
    public String getSource() { return source; }
}
