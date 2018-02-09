package uk.ac.ox.well.cortexjdk.commands.call.call;

import htsjdk.samtools.fastq.FastqRecord;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.HTreeMap;
import org.mapdb.Serializer;
import org.mapdb.serializer.SerializerArray;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.reads.Reads;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.traversal.CortexVertex;

import java.io.File;
import java.util.*;

public class AnnotateCandidates extends Module {
    @Argument(fullName="db", shortName="d", doc="Database")
    public File DB_FILE;

    @Argument(fullName="reads", shortName="r", doc="Read sequences")
    public ArrayList<File> READ_FILES;

    @Output
    public File out;

    @Override
    public void execute() {
        if (out.exists()) { out.delete(); }
        DB dbOut = DBMaker.fileDB(out).make();
        DB dbIn = DBMaker.fileDB(DB_FILE).readOnly().fileMmapEnable().make();

        HTreeMap<Integer, CortexVertex[]> dbContigsIn = initializeDbContigsInMap(dbIn);
        HTreeMap<Integer, CortexVertex[]> dbContigsOut = initializeDbContigsOutMap(dbOut);
        HTreeMap<Integer, FastqRecord[]> dbReadsEnd1Out = initializeDbReadsOutMap(dbOut, true);
        HTreeMap<Integer, FastqRecord[]> dbReadsEnd2Out = initializeDbReadsOutMap(dbOut, false);

        Map<CanonicalKmer, Set<Integer>> kmers = loadKmersInContigs(dbContigsIn, dbContigsOut);
        Map<Integer, Set<Pair<FastqRecord, FastqRecord>>> reads = annotateKmersWithReads(kmers);

        writeReadsToDatabase(reads, dbReadsEnd1Out, dbReadsEnd2Out);

        dbOut.commit();
        dbOut.close();
        dbIn.close();
    }

    private HTreeMap<Integer, FastqRecord[]> initializeDbReadsOutMap(DB dbOut, boolean firstEnd) {
        HTreeMap<Integer, FastqRecord[]> dbReads = dbOut
                .hashMap("readsEnd" + (firstEnd ? "1" : "2"))
                .keySerializer(Serializer.INTEGER)
                .valueSerializer(new SerializerArray(Serializer.JAVA))
                .counterEnable()
                .create();

        return dbReads;
    }

    @NotNull
    private HTreeMap<Integer, CortexVertex[]> initializeDbContigsOutMap(DB dbOut) {
        return (HTreeMap<Integer, CortexVertex[]>) dbOut
                    .hashMap("contigs")
                    .keySerializer(Serializer.INTEGER)
                    .valueSerializer(new SerializerArray(Serializer.JAVA))
                    .counterEnable()
                    .create();
    }

    @NotNull
    private HTreeMap<Integer, CortexVertex[]> initializeDbContigsInMap(DB dbIn) {
        return (HTreeMap<Integer, CortexVertex[]>) dbIn
                    .hashMap("contigs")
                    .keySerializer(Serializer.INTEGER)
                    .valueSerializer(new SerializerArray(Serializer.JAVA))
                    .counterEnable()
                    .open();
    }

    private Map<Integer, Set<Pair<FastqRecord, FastqRecord>>> annotateKmersWithReads(Map<CanonicalKmer, Set<Integer>> kmers) {
        Map<Integer, Set<Pair<FastqRecord, FastqRecord>>> reads = new HashMap<>();

        int kmerSize = kmers.keySet().iterator().next().length();
        int numReadsStored = 0;

        for (File readFile : READ_FILES) {
            ProgressMeter pm = new ProgressMeterFactory()
                    .header("Processing reads from " + readFile.getName() + "...")
                    .message("reads processed")
                    .updateRecord(1000000)
                    .make(log);

            Reads readsSource = new Reads(readFile);
            for (Pair<FastqRecord, FastqRecord> p : readsSource) {
                for (FastqRecord fr : Arrays.asList(p.getFirst(), p.getSecond())) {
                    if (fr != null) {
                        for (int i = 0; i <= fr.length() - kmerSize; i++) {
                            CanonicalKmer ck = new CanonicalKmer(fr.getReadString().substring(i, i + kmerSize));

                            if (kmers.containsKey(ck)) {
                                for (Integer contigIndex : kmers.get(ck)) {
                                    if (!reads.containsKey(contigIndex)) {
                                        reads.put(contigIndex, new HashSet<>());
                                    }

                                    reads.get(contigIndex).add(p);
                                    numReadsStored++;
                                }

                                break;
                            }
                        }
                    }
                }

                pm.update("reads processed (" + numReadsStored + " stored)");
            }
        }

        return reads;
    }

    @NotNull
    private Map<CanonicalKmer, Set<Integer>> loadKmersInContigs(HTreeMap<Integer, CortexVertex[]> dbContigs, HTreeMap<Integer, CortexVertex[]> dbContigsOut) {
        Map<CanonicalKmer, Set<Integer>> kmers = new HashMap<>();
        int numKmers = 0;

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Loading contigs from database...")
                .message("contigs loaded")
                .maxRecord(dbContigs.size())
                .make(log);

        for (Object o : dbContigs.keySet()) {
            int contigIndex = (Integer) o;
            List<CortexVertex> w = getWalk(dbContigs, o);

            for (CortexVertex v : w) {
                if (!kmers.containsKey(v.getCanonicalKmer())) {
                    kmers.put(v.getCanonicalKmer(), new HashSet<>());
                }

                kmers.get(v.getCanonicalKmer()).add(contigIndex);
                numKmers++;
            }

            dbContigsOut.put(contigIndex, w.toArray(new CortexVertex[w.size()]));

            pm.update();
        }

        log.info("  {} kmers loaded, {} unique", numKmers, kmers.size());

        return kmers;
    }

    private void writeReadsToDatabase(Map<Integer, Set<Pair<FastqRecord, FastqRecord>>> reads, HTreeMap<Integer, FastqRecord[]> dbReadsEnd1Out , HTreeMap<Integer, FastqRecord[]> dbReadsEnd2Out) {
        int numReads = 0;
        for (int contigIndex : reads.keySet()) {
            numReads += reads.get(contigIndex).size();
        }

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Writing records to database...")
                .message("records written")
                .maxRecord(numReads)
                .make(log);

        for (int contigIndex : reads.keySet()) {
            Set<Pair<FastqRecord, FastqRecord>> rps = reads.get(contigIndex);

            FastqRecord[] e1 = new FastqRecord[rps.size()];
            FastqRecord[] e2 = new FastqRecord[rps.size()];

            int i = 0;
            for (Pair<FastqRecord, FastqRecord> rp : rps) {
                e1[i] = rp.getFirst();
                e2[i] = rp.getSecond();

                i++;

                pm.update();
            }

            dbReadsEnd1Out.put(contigIndex, e1);
            dbReadsEnd2Out.put(contigIndex, e2);

            pm.update();
        }
    }

    private List<CortexVertex> getWalk(HTreeMap<Integer, CortexVertex[]> dbContigs, Object obIndex) {
        int contigIndex = (Integer) obIndex;

        Object[] oa = dbContigs.get(contigIndex);
        List<CortexVertex> w = new ArrayList<>();

        for (Object o : oa) {
            w.add((CortexVertex) o);
        }

        return w;
    }
}
