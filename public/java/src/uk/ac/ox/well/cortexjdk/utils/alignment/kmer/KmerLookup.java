package uk.ac.ox.well.cortexjdk.utils.alignment.kmer;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.apache.tools.ant.taskdefs.Exec;
import org.mapdb.BTreeMap;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.Serializer;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexBinaryKmer;
import uk.ac.ox.well.cortexjdk.utils.io.utils.BinaryFile;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;
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

    public static File createIndex(File reference, final int kmerSize, final String source, int nThreads, Logger log) {
        try {
            FastaSequenceFile fsf = new FastaSequenceFile(reference, true);
            ReferenceSequence rseq;

            ExecutorService executor = Executors.newFixedThreadPool(nThreads);
            //List<Callable<Map<Integer, Map<CortexBinaryKmer, List<Integer>>>>> results = new ArrayList<>();
            List<Future<Map<Integer, Map<CortexBinaryKmer, List<Integer>>>>> results = new ArrayList<>();

            while ((rseq = fsf.nextSequence()) != null) {
                final int contigIndex = rseq.getContigIndex();
                final String rname = rseq.getName();
                final byte[] seq = rseq.getBases();

                Callable<Map<Integer, Map<CortexBinaryKmer, List<Integer>>>> task = () -> {
                    String threadName = Thread.currentThread().getName();

                    ProgressMeter pm = new ProgressMeterFactory()
                            .header("Processing contig (" + threadName + ", " + rname + ")")
                            .message("processed (" + threadName + ", " + rname + ")")
                            .maxRecord(seq.length - kmerSize)
                            .updateRecord((seq.length - kmerSize)/2)
                            .make(log);

                    Map<Integer, Map<CortexBinaryKmer, List<Integer>>> kmerPos = new TreeMap<>();
                    kmerPos.put(contigIndex, new TreeMap<>());

                    for (int i = 0; i <= seq.length - kmerSize; i++) {
                        byte[] kmer = new byte[kmerSize];
                        System.arraycopy(seq, i, kmer, 0, kmerSize);

                        if (SequenceUtils.isValidNucleotideSequence(kmer)) {
                            CortexBinaryKmer cbk = new CortexBinaryKmer(kmer);

                            if (!kmerPos.get(contigIndex).containsKey(cbk)) {
                                kmerPos.get(contigIndex).put(cbk, new ArrayList<>());
                            }

                            kmerPos.get(contigIndex).get(cbk).add(i);
                        }

                        pm.update();
                    }

                    return kmerPos;
                };

                results.add(executor.submit(task));
            }

            try {
                Map<CortexBinaryKmer, List<int[]>> kmers = new TreeMap<>();
                boolean isFinished;

                do {
                    isFinished = true;
                    for (Future<Map<Integer, Map<CortexBinaryKmer, List<Integer>>>> future : results) {
                        if (future.isDone()) {
                            for (int contigIndex : future.get().keySet()) {
                                for (CortexBinaryKmer cbk : future.get().get(contigIndex).keySet()) {
                                    if (!kmers.containsKey(cbk)) {
                                        kmers.put(cbk, new ArrayList<>());
                                    }

                                    for (int p : future.get().get(contigIndex).get(cbk)) {
                                        int[] l = new int[]{ contigIndex, p };

                                        kmers.get(cbk).add(l);
                                    }
                                }
                            }
                        }

                        isFinished &= future.isDone();
                    }
                } while (!isFinished);

                executor.shutdownNow();

                File dbFile = new File(reference.getAbsolutePath() + ".idx");
                BinaryFile bf = new BinaryFile(dbFile, "rw");

                storeAtoms(bf, 1, kmerSize, source);
                storeKmers(bf, kmers, log);

                /*
                for (CortexBinaryKmer cbk : kmers.keySet()) {
                    log.info("{} {}", cbk, kmers.get(cbk));
                }
                */

                bf.close();

                return dbFile;
            } catch (InterruptedException e) {
                throw new CortexJDKException("Interrupted");
            } catch (ExecutionException e) {
                throw new CortexJDKException("Execution problem");
            }
        } catch (FileNotFoundException e) {
            throw new CortexJDKException("File not found", e);
        } catch (IOException e) {
            throw new CortexJDKException("IOException", e);
        }
    }

    private static void storeAtoms(BinaryFile bf, int version, int kmerSize, String source) throws IOException {
        bf.writeBytes("KMRIDX");
        bf.writeInt(version);
        bf.writeInt(kmerSize);
        bf.writeInt(source.length());
        bf.writeBytes(source);
        bf.writeBytes("KMRIDX");
    }

    private static void storeKmers(BinaryFile bf, Map<CortexBinaryKmer, List<int[]>> kmers, Logger log) throws IOException {
        ProgressMeter pm = new ProgressMeterFactory()
                .header("Writing index")
                .message("records written")
                .maxRecord(kmers.size())
                .make(log);

        int numLoci = 0;
        for (CortexBinaryKmer cbk : kmers.keySet()) {
            for (int[] l : kmers.get(cbk)) {
                long[] k = cbk.getBinaryKmer();

                for (int i = 0; i < k.length; i++) {
                    bf.writeLong(k[i]);
                }

                bf.writeInt(l[0]);
                bf.writeInt(l[1]);

                numLoci++;
            }

            pm.update();
        }

        log.info("Wrote {} kmers from {} loci", kmers.size(), numLoci);
    }
}
