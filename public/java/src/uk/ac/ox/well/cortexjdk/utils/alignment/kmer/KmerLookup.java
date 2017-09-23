package uk.ac.ox.well.cortexjdk.utils.alignment.kmer;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexBinaryKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.utils.BinaryFile;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;

import org.slf4j.Logger;

public class KmerLookup {
    private File refFile;
    private IndexedFastaSequenceFile ref;
    private KmerStream idx;

    public KmerLookup(File refFile) { initialize(refFile); }

    private void initialize(File refFile) {
        try {
            this.refFile = refFile;
            this.ref = new IndexedFastaSequenceFile(refFile);

            if (this.ref.getSequenceDictionary() == null) {
                throw new CortexJDKException("Reference must have a sequence dictionary");
            }

            idx = new KmerStream(new File(refFile.getAbsolutePath() + ".kidx"));
        } catch (FileNotFoundException e) {
            throw new CortexJDKException("File not found", e);
        }
    }

    public int getKmerSize() { return idx.getKmerSize(); }

    public String getSource() { return idx.getSource(); }

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

    public Set<Interval> findKmer(String seq) {
        if (seq.length() < idx.getKmerSize()) { throw new CortexJDKException("Minimum kmer size for lookups is " + idx.getKmerSize()); }

        List<List<Interval>> allIntervals = new ArrayList<>();
        for (int i = 0; i <= seq.length() - idx.getKmerSize(); i++) {
            String sk = seq.substring(i, i + idx.getKmerSize());
            List<KmerStreamRec> recs = idx.find(sk);

            Set<Interval> intervals = new HashSet<>();
            for (KmerStreamRec ksr : recs) {
                String contigName = ref.getSequenceDictionary().getSequence(ksr.contigIndex).getSequenceName();

                String refk = ref.getSubsequenceAt(contigName, ksr.offset + 1, ksr.offset + idx.getKmerSize()).getBaseString();
                boolean isNegative = refk.equals(SequenceUtils.reverseComplement(sk));

                intervals.add(new Interval(contigName, ksr.offset + 1, ksr.offset + idx.getKmerSize(), isNegative, null));
            }

            allIntervals.add(new ArrayList<>(intervals));
        }

        return combineIntervals(allIntervals, seq);
    }

    private Set<Interval> combineIntervals(List<List<Interval>> allIntervals, String seq) {
        List<Interval> lastIntervals = allIntervals.get(0);

        for (int i = 1; i < allIntervals.size(); i++) {
            List<Interval> newIntervals = new ArrayList<>();

            for (int j = 0; j < allIntervals.get(i).size(); j++) {
                Interval thisInterval = allIntervals.get(i).get(j);

                for (Interval lastInterval : lastIntervals) {
                    if (lastInterval.getIntersectionLength(thisInterval) == thisInterval.length() - 1 &&
                        lastInterval.isNegativeStrand() == thisInterval.isNegativeStrand()) {

                        Interval combinedInterval = new Interval(
                                lastInterval.getContig(),
                                lastInterval.getStart() < thisInterval.getStart() ? lastInterval.getStart() : thisInterval.getStart(),
                                lastInterval.getEnd() > thisInterval.getEnd() ? lastInterval.getEnd() : thisInterval.getEnd(),
                                lastInterval.isNegativeStrand(),
                                null
                        );

                        newIntervals.add(combinedInterval);
                    }
                }
            }

            lastIntervals = newIntervals;
        }

        Set<Interval> finalIntervals = new HashSet<>();
        for (Interval interval : lastIntervals) {
            if (interval.length() - getKmerSize() + 1 == allIntervals.size()) {
                String fwd = ref.getSubsequenceAt(interval.getContig(), interval.getStart(), interval.getEnd()).getBaseString();
                String rev = SequenceUtils.reverseComplement(fwd);

                String nseq = interval.isPositiveStrand() ? fwd : rev;

                if (nseq.equals(seq)) {
                    finalIntervals.add(interval);
                }
            }
        }

        return finalIntervals;
    }

    public static File createIndex(File reference, final int kmerSize, final String source, int nThreads, Logger log) {
        try {
            List<KmerStream> files = getSortedKmerStreams(reference, kmerSize, nThreads, log);

            File dbFile = new File(reference.getAbsolutePath() + ".kidx");
            BinaryFile bf = new BinaryFile(dbFile, "rw");

            storeAtoms(bf, kmerSize, source);
            storeKmers(bf, files, log);

            bf.close();

            return dbFile;
        } catch (FileNotFoundException e) {
            throw new CortexJDKException("File not found", e);
        } catch (IOException e) {
            throw new CortexJDKException("IOException", e);
        }
    }

    @NotNull
    private static List<KmerStream> getSortedKmerStreams(final File reference, int kmerSize, int nThreads, Logger log) {
        FastaSequenceFile fsf = new FastaSequenceFile(reference, true);
        ReferenceSequence rseq;

        ExecutorService executor = Executors.newFixedThreadPool(nThreads);
        List<Future<Pair<File, Integer>>> results = new ArrayList<>();

        while ((rseq = fsf.nextSequence()) != null) {
            final int contigIndex = rseq.getContigIndex();
            final String rname = rseq.getName();
            final byte[] seq = rseq.getBases();

            Callable<Pair<File, Integer>> task = () -> {
                String threadName = Thread.currentThread().getName();

                ProgressMeter pm = new ProgressMeterFactory()
                    .header("Processing contig (" + threadName + ", " + rname + ")")
                    .message("processed (" + threadName + ", " + rname + ")")
                    .maxRecord(seq.length - kmerSize)
                    .make(log);

                List<Pair<CortexBinaryKmer, Integer>> kmerPos = new ArrayList<>(seq.length - kmerSize);

                for (int i = 0; i <= seq.length - kmerSize; i++) {
                    byte[] kmer = new byte[kmerSize];
                    System.arraycopy(seq, i, kmer, 0, kmerSize);

                    if (SequenceUtils.isValidNucleotideSequence(kmer)) {
                        CortexBinaryKmer cbk = new CortexBinaryKmer(kmer);

                        kmerPos.add(new Pair<>(cbk, i));
                    }

                    pm.update();
                }

                if (log != null) { log.info("  sorting kmers ({}, {})", threadName, rname); }

                kmerPos.sort((o1, o2) -> {
                    if (o1.getFirst().equals(o2.getFirst())) {
                        return o1.getSecond().compareTo(o2.getSecond());
                    }

                    return o1.getFirst().compareTo(o2.getFirst());
                });

                File f = File.createTempFile(reference.getName(), ".ks", reference.getAbsoluteFile().getParentFile());
                f.deleteOnExit();

                if (log != null) { log.info("  writing kmers to {} ({}, {})", f.getAbsolutePath(), threadName, rname); }

                BinaryFile bf = new BinaryFile(f, "rw");

                storeAtoms(bf, kmerSize, rname);
                storeKmers(bf, kmerPos, contigIndex);

                bf.close();

                return new Pair<>(f, contigIndex);
            };

            results.add(executor.submit(task));
        }

        boolean isFinished;
        do {
            isFinished = true;
            for (Future<Pair<File, Integer>> future : results) {
                isFinished &= future.isDone();
            }
        } while (!isFinished);

        executor.shutdownNow();

        List<KmerStream> files = new ArrayList<>();
        results.forEach(f -> {
            try {
                files.add(new KmerStream(f.get().getFirst()));
            } catch (InterruptedException e) {
                throw new CortexJDKException("Interrupted", e);
            } catch (ExecutionException e) {
                throw new CortexJDKException("Execution problem", e);
            }
        });

        return files;
    }

    private static void storeAtoms(BinaryFile bf, int kmerSize, String source) throws IOException {
        bf.writeBytes("KMRIDX");
        bf.writeInt(1); // version
        bf.writeInt(kmerSize);
        bf.writeInt(source.length());
        bf.writeBytes(source);
        bf.writeBytes("KMRIDX");
    }

    private static void storeKmers(BinaryFile bf, List<KmerStream> files, Logger log) throws IOException {
        ProgressMeter pm = new ProgressMeterFactory()
                .header("Writing index")
                .message("records written")
                .updateRecord(files.get(0).getNumRecords() / 10)
                .make(log);

        KmerStreamRec lowestKmer;
        do {
            lowestKmer = null;
            for (int i = 0; i < files.size(); i++) {
                KmerStreamRec ksr = files.get(i).peek();

                if (lowestKmer == null || ksr.kmerCompare(lowestKmer) < 0) {
                    lowestKmer = ksr;
                }
            }

            if (lowestKmer != null) {
                for (int i = 0; i < files.size(); i++) {
                    KmerStreamRec ksr = files.get(i).peek();

                    while (ksr != null && ksr.kmerCompare(lowestKmer) == 0) {
                        long[] l = ksr.cbk.getBinaryKmer();
                        for (int j = 0; j < l.length; j++) {
                            bf.writeLong(l[j]);
                        }
                        bf.writeInt(ksr.contigIndex);
                        bf.writeInt(ksr.offset);

                        files.get(i).advance();
                        ksr = files.get(i).peek();

                        pm.update();
                    }
                }
            }
        } while (lowestKmer != null);

        if (log != null) { log.info("Wrote {} records", pm.pos()); }
    }

    private static void storeKmers(BinaryFile bf, List<Pair<CortexBinaryKmer, Integer>> kmerPos, int contigIndex) {
        try {
            for (Pair<CortexBinaryKmer, Integer> p : kmerPos) {
                for (long l : p.getFirst().getBinaryKmer()) {
                    bf.writeLong(l);
                }
                bf.writeInt(contigIndex);
                bf.writeInt(p.getSecond());
            }
        } catch (IOException e) {
            throw new CortexJDKException("IOException", e);
        }
    }
}
