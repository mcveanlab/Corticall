package uk.ac.ox.well.cortexjdk.utils.io.cortex.links;

import org.jetbrains.annotations.NotNull;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.HTreeMap;
import org.mapdb.Serializer;
import org.mapdb.serializer.SerializerArrayTuple;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;

import java.io.File;
import java.util.*;

public class CortexLinks implements Map<CortexKmer, CortexLinksRecord> {
    private CortexLinksMap clm = null;

    private DB db = null;
    private NavigableSet<Object[]> linkIndex;
    private HTreeMap<Integer, String> linkSources;
    private HTreeMap<Integer, String> linkColors;
    private Map<String, Integer> linkSamples;

    public CortexLinks(String linksFilePath) {
        loadLinks(new File(linksFilePath));
    }

    public CortexLinks(File linksFile) {
        loadLinks(linksFile);
    }

    private void loadLinks(File linksFile) {
        File dbFile = linksFile.getAbsolutePath().endsWith(".linkdb") ? linksFile : new File(linksFile.getAbsolutePath() + ".linkdb");

        if (dbFile.exists()) {
            db = DBMaker
                    .fileDB(dbFile)
                    .fileMmapEnable()
                    .readOnly()
                    .make();

            int version = db.atomicInteger("version").open().get();
            if (version != 1) {
                throw new CortexJDKException("Expected linkdb file of version=1, found version=" + version);
            }

            linkColors = db.hashMap("colors")
                    .keySerializer(Serializer.INTEGER)
                    .valueSerializer(Serializer.STRING)
                    .counterEnable()
                    .open();

            linkSources = db.hashMap("sources")
                    .keySerializer(Serializer.INTEGER)
                    .valueSerializer(Serializer.STRING)
                    .counterEnable()
                    .open();

            linkSamples = new HashMap<>();
            for (Object o : linkColors.keySet()) {
                Integer c = (Integer) o;
                String sampleName = linkColors.get(c);
                linkSamples.put(sampleName, c);
            }

            linkIndex = db.treeSet("links")
                    .serializer(new SerializerArrayTuple(Serializer.BYTE_ARRAY, Serializer.BYTE_ARRAY, Serializer.BYTE_ARRAY))
                    .counterEnable()
                    .open();
        } else {
            clm = new CortexLinksMap(linksFile);
        }
    }

    public int getNumSources() {
        return linkSources == null ? 0 : linkSources.size();
    }

    public String getSourceForIndex(int index) {
        if (linkSources == null) {
            return "unknown";
        }

        return linkSources.get(index);
    }

    public Collection<String> getSources() {
        if (linkSources == null) {
            return new ArrayList<>();
        }

        return linkSources.values();
    }

    public Integer getColorForSampleName(String sampleName) {
        if (linkSamples == null) {
            return clm.getCortexGraphLinks().getColorForSampleName(sampleName);
        }

        return linkSamples.get(sampleName);
    }

    public String getSampleNameForColor(int c) {
        if (linkColors == null) {
            return clm.getCortexGraphLinks().getColor(c).getSampleName();
        }

        return linkColors.get(c);
    }

    @Override
    public int size() {
        return clm == null ? linkIndex.size() : clm.size();
    }

    @Override
    public boolean isEmpty() {
        return clm == null ? linkIndex.size() == 0 : clm.size() == 0;
    }

    public boolean containsKey(byte[] bk) {
        if (clm == null) {
            Set subset = linkIndex.subSet(
                    new Object[]{bk},
                    new Object[]{bk, null, null});

            return subset.size() > 0;
        } else {
            return clm.containsKey(bk);
        }
    }

    public boolean containsKey(String sk) {
        return containsKey(sk.getBytes());
    }

    public boolean containsKey(CortexKmer ck) {
        return containsKey(ck.getKmerAsBytes());
    }

    public CortexLinksRecord get(byte[] bk) {
        if (clm == null) {
            Set subset = linkIndex.subSet(
                    new Object[]{bk},
                    new Object[]{bk, null, null});

            List<CortexJunctionsRecord> cjrs = new ArrayList<>();
            for (Object o : subset) {
                Object[] oa = (Object[]) o;
                CortexJunctionsRecord cjr = decode((byte[]) oa[1]);
                cjrs.add(cjr);
            }

            return new CortexLinksRecord(new String(bk), cjrs);
        } else {
            return clm.get(new String(bk));
        }
    }

    public CortexLinksRecord get(String sk) { return get(sk.getBytes()); }

    public CortexLinksRecord get(CortexKmer ck) { return get(ck.getKmerAsBytes()); }

    @Override
    public CortexLinksRecord get(Object key) {
        throw new CortexJDKException("Not implemented.");
    }

    @Override
    public boolean containsKey(Object key) {
        throw new CortexJDKException("Not implemented.");
    }

    @Override
    public boolean containsValue(Object value) {
        throw new CortexJDKException("Not implemented.");
    }

    @Override
    public CortexLinksRecord put(CortexKmer key, CortexLinksRecord value) { throw new CortexJDKException("Not implemented."); }

    @Override
    public CortexLinksRecord remove(Object key) {
        throw new CortexJDKException("Not implemented.");
    }

    @Override
    public void putAll(@NotNull Map<? extends CortexKmer, ? extends CortexLinksRecord> m) { throw new CortexJDKException("Not implemented."); }

    @Override
    public void clear() {
        throw new CortexJDKException("Not implemented.");
    }

    @NotNull
    @Override
    public Set<CortexKmer> keySet() {
        if (clm == null) {
            Set<CortexKmer> keys = new HashSet<>();

            for (Object[] o : linkIndex) {
                byte[] o0 = (byte[]) o[0];
                CortexKmer ck = new CortexKmer(o0);

                keys.add(ck);
            }

            return keys;
        } else {
            return clm.keySet();
        }
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
        throw new CortexJDKException("Not implemented.");
    }

    @NotNull
    @Override
    public Set<Entry<CortexKmer, CortexLinksRecord>> entrySet() {
        throw new CortexJDKException("Not implemented.");
    }
}
