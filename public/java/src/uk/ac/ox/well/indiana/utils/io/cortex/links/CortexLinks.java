package uk.ac.ox.well.indiana.utils.io.cortex.links;

import org.jetbrains.annotations.NotNull;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.HTreeMap;
import org.mapdb.Serializer;
import org.mapdb.serializer.SerializerArrayTuple;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;

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
                throw new IndianaException("Expected linkdb file of version=1, found version=" + version);
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
        return linkSources.size();
    }

    public String getSourceForIndex(int index) {
        return linkSources.get(index);
    }

    public Collection<String> getSources() {
        return linkSources.values();
    }

    public Integer getColorForSampleName(String sampleName) { return linkSamples.get(sampleName); }

    public String getSampleNameForColor(int c) { return linkColors.get(c); }

    @Override
    public int size() {
        return clm == null ? linkIndex.size() : clm.size();
    }

    @Override
    public boolean isEmpty() {
        return clm == null ? linkIndex.size() == 0 : clm.size() == 0;
    }

    public boolean containsKey(String sk) {
        if (clm == null) {
            Set subset = linkIndex.subSet(
                    new Object[]{sk.getBytes()},
                    new Object[]{sk.getBytes(), null, null});

            return subset.size() > 0;
        } else {
            return clm.containsKey(sk);
        }
    }

    public boolean containsKey(CortexKmer ck) {
        if (clm == null) {
            Set subset = linkIndex.subSet(
                    new Object[]{ck.getKmerAsBytes()},
                    new Object[]{ck.getKmerAsBytes(), null, null});

            return subset.size() > 0;
        } else {
            return clm.containsKey(ck);
        }
    }

    public CortexLinksRecord get(String sk) {
        if (clm == null) {
            Set subset = linkIndex.subSet(
                    new Object[]{sk.getBytes()},
                    new Object[]{sk.getBytes(), null, null});


            List<CortexJunctionsRecord> cjrs = new ArrayList<>();
            for (Object o : subset) {
                Object[] oa = (Object[]) o;
                CortexJunctionsRecord cjr = decode((byte[]) oa[1]);
                cjrs.add(cjr);
            }

            return new CortexLinksRecord(sk, cjrs);
        } else {
            return clm.get(sk);
        }
    }

    public CortexLinksRecord get(CortexKmer ck) {
        if (clm == null) {
            Set subset = linkIndex.subSet(
                    new Object[]{ck.getKmerAsBytes()},
                    new Object[]{ck.getKmerAsBytes(), null, null});


            List<CortexJunctionsRecord> cjrs = new ArrayList<>();
            for (Object o : subset) {
                Object[] oa = (Object[]) o;
                CortexJunctionsRecord cjr = decode((byte[]) oa[1]);
                cjrs.add(cjr);
            }

            return new CortexLinksRecord(ck.getKmerAsString(), cjrs);
        } else {
            return clm.get(ck);
        }
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
        throw new IndianaException("Not implemented.");
    }

    @NotNull
    @Override
    public Set<Entry<CortexKmer, CortexLinksRecord>> entrySet() {
        throw new IndianaException("Not implemented.");
    }
}
