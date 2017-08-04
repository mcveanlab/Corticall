package uk.ac.ox.well.indiana.utils.io.cortex.links;

import org.jetbrains.annotations.NotNull;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.Serializer;
import org.mapdb.serializer.SerializerArrayTuple;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexBinaryKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.File;
import java.util.*;

public class CortexLinks implements Map<CortexKmer, CortexLinksRecord> {
    private File dbFile;
    private DB db;
    private NavigableSet<Object[]> linkIndex;

    public CortexLinks(String cortexLinksFilePath) {
        loadLinksDB(new File(cortexLinksFilePath + ".linkdb"));
    }

    public CortexLinks(File cortexLinksFile) {
        loadLinksDB(new File(cortexLinksFile.getAbsolutePath() + ".linkdb"));
    }

    private void loadLinksDB(File dbFile) {
        this.dbFile = dbFile;

        db = DBMaker
                .fileDB(this.dbFile)
                .fileMmapEnable()
                .readOnly()
                .make();

        linkIndex = db.treeSet("links")
                .serializer(new SerializerArrayTuple(Serializer.BYTE_ARRAY, Serializer.BYTE_ARRAY))
                .counterEnable()
                .counterEnable()
                .counterEnable()
                .open();
    }

    @Override
    public int size() {
        return linkIndex.size();
    }

    @Override
    public boolean isEmpty() {
        return linkIndex.size() == 0;
    }

    public boolean containsKey(String sk) {
        Set subset = linkIndex.subSet(
                new Object[]{sk.getBytes()},
                new Object[]{sk.getBytes(), null});

        return subset.size() > 0;
    }

    public boolean containsKey(CortexKmer ck) {
        Set subset = linkIndex.subSet(
                new Object[]{ck.getKmerAsBytes()},
                new Object[]{ck.getKmerAsBytes(), null});

        return subset.size() > 0;
    }

    public CortexLinksRecord get(String sk) {
        Set subset = linkIndex.subSet(
                new Object[]{sk.getBytes()},
                new Object[]{sk.getBytes(), null});


        List<CortexJunctionsRecord> cjrs = new ArrayList<>();
        for (Object o : subset) {
            Object[] oa = (Object[]) o;
            CortexJunctionsRecord cjr = decode((byte[]) oa[1]);
            cjrs.add(cjr);
        }

        return new CortexLinksRecord(sk, cjrs);
    }

    public CortexLinksRecord get(CortexKmer ck) {
        Set subset = linkIndex.subSet(
                new Object[]{ck.getKmerAsBytes()},
                new Object[]{ck.getKmerAsBytes(), null});


        List<CortexJunctionsRecord> cjrs = new ArrayList<>();
        for (Object o : subset) {
            Object[] oa = (Object[]) o;
            CortexJunctionsRecord cjr = decode((byte[]) oa[1]);
            cjrs.add(cjr);
        }

        return new CortexLinksRecord(ck.getKmerAsString(), cjrs);
    }

    @Override
    public CortexLinksRecord get(Object key) {
        throw new IndianaException("Not implemented.");
    }

    @Override
    public boolean containsKey(Object key) {
        throw new IndianaException("Not implemented.");
    }

    @Override
    public boolean containsValue(Object value) {
        throw new IndianaException("Not implemented.");
    }

    @Override
    public CortexLinksRecord put(CortexKmer key, CortexLinksRecord value) {
        throw new IndianaException("Not implemented.");
    }

    @Override
    public CortexLinksRecord remove(Object key) {
        throw new IndianaException("Not implemented.");
    }

    @Override
    public void putAll(@NotNull Map<? extends CortexKmer, ? extends CortexLinksRecord> m) {
        throw new IndianaException("Not implemented.");
    }

    @Override
    public void clear() {
        throw new IndianaException("Not implemented.");
    }

    @NotNull
    @Override
    public Set<CortexKmer> keySet() {
        Set<CortexKmer> keys = new HashSet<>();

        for (Object[] o : linkIndex) {
            byte[] o0 = (byte[]) o[0];
            CortexKmer ck = new CortexKmer(o0);

            keys.add(ck);
        }

        return keys;
    }

    private CortexJunctionsRecord decode(byte[] bline) {
        String line = new String(bline);
        String[] linkLine = line.split("\\s+");

        String orientation = linkLine[0];
        int numJunctions = Integer.valueOf(linkLine[1]);

        String[] covs = linkLine[2].split(",");
        int[] coverages = new int[covs.length];
        for (int i = 0; i < covs.length; i++) {
            coverages[i] = Integer.valueOf(covs[i]);
        }

        String junctions = linkLine[3];

        return new CortexJunctionsRecord(orientation.equals("F"), -1, numJunctions, coverages, junctions);
    }

    @NotNull
    @Override
    public Collection<CortexLinksRecord> values() {
        throw new IndianaException("Not implemented.");
    }

    @NotNull
    @Override
    public Set<Entry<CortexKmer, CortexLinksRecord>> entrySet() {
        throw new IndianaException("Not implemented.");
    }
}
