package uk.ac.ox.well.cortexjdk.utils.alignment.kmer;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.mapdb.BTreeMap;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.Serializer;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexBinaryKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.slf4j.Logger;

public class KmerLookup {
    private File refFile;
    private IndexedFastaSequenceFile ref;
    private String source;
    private Map<Integer, BTreeMap<long[], int[]>> kmerIndices = new HashMap<>();

    public KmerLookup(File refFile) { initialize(refFile); }

    private void initialize(File refFile) {
        try {
            this.refFile = refFile;
            this.ref = new IndexedFastaSequenceFile(refFile);

            if (this.ref.getSequenceDictionary() == null) {
                throw new CortexJDKException("Reference must have a sequence dictionary");
            }

            loadIndices(refFile);
        } catch (FileNotFoundException e) {
            throw new CortexJDKException("File not found", e);
        }
    }

    private void loadIndices(File refFile) {
        File[] matchingFiles = refFile.getParentFile().listFiles((d, n) -> n.startsWith(refFile.getName()) && n.endsWith(".kmerdb"));

        if (matchingFiles == null || matchingFiles.length == 0) {
            throw new CortexJDKException("No .kmerdb files found for fasta '" + refFile.getAbsolutePath() + "'");
        }

        Pattern p = Pattern.compile(".k(\\d+).kmerdb");

        for (File dbFile : matchingFiles) {
            Matcher m = p.matcher(dbFile.getName());

            if (m.find()) {
                int kmerSize = Integer.valueOf(m.group(1));

                DB db = DBMaker.fileDB(dbFile)
                        .fileMmapEnable()
                        .closeOnJvmShutdown()
                        .readOnly()
                        .make();

                int version = db.atomicInteger("version").open().get();
                if (version != 1) {
                    throw new CortexJDKException("Expected .kmerdb file of version=1, found version=" + version);
                }

                String source = db.atomicString("source").open().get();
                if (this.source != null && !this.source.equals(source)) {
                    throw new CortexJDKException("Source specified in '" + refFile.getAbsolutePath() + "' does not match previous .kmerdb sources.");
                }
                this.source = source;

                kmerIndices.put(kmerSize, db.treeMap("index")
                        .keySerializer(Serializer.LONG_ARRAY)
                        .valueSerializer(Serializer.INT_ARRAY)
                        .counterEnable()
                        .open());
            }
        }
    }

    public Set<Integer> getKmerSizes() { return kmerIndices.keySet(); }

    public String getSource() { return source; }

    public IndexedFastaSequenceFile getReferenceSequence() {
        return ref;
    }

    public String findKmer(Interval interval) {
        if (ref.getSequenceDictionary().getSequenceIndex(interval.getContig()) == -1) {
            throw new CortexJDKException("Contig '" + interval.getContig() + "' was not found in reference '" + refFile.getAbsolutePath() + "'");
        }

        if (interval.getStart() > 0 && interval.getEnd() <= ref.getSequenceDictionary().getSequence(interval.getContig()).getSequenceLength()) {
            ReferenceSequence rseq = ref.getSubsequenceAt(interval.getContig(), interval.getStart(), interval.getEnd());

            if (rseq != null) {
                if (interval.isPositiveStrand()) {
                    return rseq.getBaseString();
                } else {
                    return SequenceUtils.reverseComplement(rseq.getBaseString());
                }
            }
        }

        return null;
    }

    public Set<Interval> findKmer(String sk) {
        if (!kmerIndices.containsKey(sk.length())) {
            throw new CortexJDKException("No index at k=" + sk.length() + " for '" + refFile.getAbsolutePath() + "'");
        }

        Set<Interval> intervals = new HashSet<>();

        int[] l = kmerIndices.get(sk.length()).get(new CortexBinaryKmer(sk.getBytes()).getBinaryKmer());

        if (l != null) {
            for (int i = 0; i < l.length - 1; i += 2) {
                String chr = ref.getSequenceDictionary().getSequence(l[i]).getSequenceName();
                int pos = l[i + 1];

                String fw = ref.getSubsequenceAt(chr, pos + 1, pos + sk.length()).getBaseString();
                String rc = SequenceUtils.reverseComplement(fw);

                if (sk.equals(fw)) {
                    intervals.add(new Interval(chr, pos + 1, pos + sk.length()));
                } else if (sk.equals(rc)) {
                    intervals.add(new Interval(chr, pos + 1, pos + sk.length(), true, null));
                }
            }
        }

        return intervals;
    }

    public static File createIndex(File reference, String source, int kmerSize, Logger log) {
        File dbFile = new File(reference.getAbsoluteFile() + ".k" + kmerSize + ".kmerdb");

        DB db = DBMaker
                .fileDB(dbFile)
                .fileMmapEnable()
                .make();

        storeAtoms(db, 1, source, kmerSize);
        storeKmers(db, reference, kmerSize, log);

        db.close();

        return dbFile;
    }

    private static void storeAtoms(DB db, int version, String source, int kmerSize) {
        db.atomicInteger("version", version).create();
        db.atomicInteger("kmerSize", kmerSize).create();
        db.atomicString("source", source).create();

        db.commit();
    }

    private static void storeKmers(DB db, File reference, int kmerSize, Logger log) {
        BTreeMap<long[], int[]> kmerIndex = db.treeMap("index")
                .keySerializer(Serializer.LONG_ARRAY)
                .valueSerializer(Serializer.INT_ARRAY)
                .create();

        FastaSequenceFile ref = new FastaSequenceFile(reference, true);
        ReferenceSequence rseq;
        while ((rseq = ref.nextSequence()) != null) {
            if (log != null) { log.info("  {}", rseq.getName()); }

            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                String sk = seq.substring(i, i + kmerSize);

                if (sk.matches("^[ACGTacgt]+$")) {
                    CortexBinaryKmer cbk = new CortexBinaryKmer(sk.getBytes());

                    int[] l = {rseq.getContigIndex(), i};

                    if (!kmerIndex.containsKey(cbk.getBinaryKmer())) {
                        kmerIndex.put(cbk.getBinaryKmer(), l);
                    } else {
                        int[] l1 = kmerIndex.get(cbk.getBinaryKmer());
                        int[] l2 = new int[l1.length + 2];

                        System.arraycopy(l1, 0, l2, 0, l1.length);
                        System.arraycopy(l, 0, l2, l1.length, l.length);

                        kmerIndex.put(cbk.getBinaryKmer(), l2);
                    }
                }
            }

            db.commit();
        }
    }
}
