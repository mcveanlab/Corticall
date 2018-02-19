package uk.ac.ox.well.cortexjdk.utils.io.reads;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.utils.LineReader;

import java.io.File;
import java.util.Iterator;

import static uk.ac.ox.well.cortexjdk.utils.io.reads.Reads.SeqType.*;

public class Reads implements Iterable<Pair<FastqRecord, FastqRecord>>, Iterator<Pair<FastqRecord, FastqRecord>> {
    public enum SeqType { FASTQ, SAM, FASTA, TEXT }

    private File[] files;
    private FastqReader[] fqends;
    private SAMRecordIterator[] sams;
    private FastaSequenceFile[] faends;
    private LineReader[] lends;

    private SeqType type;
    private int numEnds = 0;

    private Pair<FastqRecord, FastqRecord> nextRecord;

    public Reads(File readsFile) {
        String[] pieces = readsFile.getAbsolutePath().split(":");
        numEnds = pieces.length;

        if (numEnds > 2) {
            throw new CortexJDKException("Too many files provided for paired-end mode (" + Joiner.on(", ").join(pieces));
        }

        files = new File[numEnds];

        for (int i = 0; i < numEnds; i++) {
            String readsPath = pieces[i];

            if (! new File(readsPath).exists()) {
                throw new CortexJDKException("File not found: '" + readsPath + "'");
            } else {
                files[i] = new File(readsPath);
            }

            if (readsPath.endsWith(".fastq.gz") || readsPath.endsWith(".fq.gz") || readsPath.endsWith(".fastq") || readsPath.endsWith(".fq")) {
                if (fqends == null) {
                    fqends = new FastqReader[numEnds];
//                    fqends = new FqReader[numEnds];
                    type = FASTQ;
                } else {
                    if (type != FASTQ) {
                        throw new CortexJDKException("Paired-end files must be of the same file type (" + pieces[i-1] + ", " + readsPath + ")");
                    }
                }

                fqends[i] = new FastqReader(new File(readsPath));
//                fqends[i] = new FqReader(new File(readsPath));
            } else if (readsPath.endsWith(".sam") || readsPath.endsWith(".bam")) {
                if (sams == null) {
                    sams = new SAMRecordIterator[numEnds];
                    type = SAM;
                } else {
                    if (type != SAM) {
                        throw new CortexJDKException("Paired-end files must be of the same file type (" + pieces[i-1] + ", " + readsPath + ")");
                    }
                }

                sams[i] = SamReaderFactory
                             .make()
                             .open(new File(readsPath))
                             .iterator();
            } else if (readsPath.endsWith(".fasta") || readsPath.endsWith(".fa") || readsPath.endsWith(".fna")) {
                if (faends == null) {
                    faends = new FastaSequenceFile[numEnds];
                    type = FASTA;
                } else {
                    if (type != FASTA) {
                        throw new CortexJDKException("Paired-end files must be of the same file type (" + pieces[i-1] + ", " + readsPath + ")");
                    }
                }

                faends[i] = new FastaSequenceFile(new File(readsPath), false);
            } else {
                if (lends == null) {
                    lends = new LineReader[numEnds];
                    type = TEXT;
                } else {
                    if (type != TEXT) {
                        throw new CortexJDKException("Paired-end files must be of the same file type (" + pieces[i-1] + ", " + readsPath + ")");
                    }
                }

                lends[i] = new LineReader(new File(readsPath));
            }
        }

        nextRecord = getNextRecord();
    }

    private Pair<FastqRecord, FastqRecord> getNextRecord() {
        FastqRecord[] ends = new FastqRecord[2];

        for (int i = 0; i < numEnds; i++) {
            switch (type) {
                case FASTQ:
                    if (fqends[i].hasNext()) {
                        ends[i] = fqends[i].next();
                    }
                    break;
                case SAM:
                    if (sams[i].hasNext()) {
                        SAMRecord sr = sams[i].next();
                        ends[i] = new FastqRecord(sr.getReadName(), sr.getReadString(), null, sr.getBaseQualityString());
                    }
                    break;
                case FASTA:
                    ReferenceSequence rseq = faends[i].nextSequence();
                    if (rseq != null) {
                        ends[i] = new FastqRecord(rseq.getName(), rseq.getBaseString(), null, StringUtil.repeatCharNTimes('I', rseq.length()));
                    }
                    break;
                case TEXT:
                    if (lends[i].hasNext()) {
                        String line = lends[i].getNextRecord();
                        ends[i] = new FastqRecord(null, line, null, StringUtil.repeatCharNTimes('I', line.length()));
                    }
                    break;
            }
        }

        return ends[0] == null ? null : new Pair<>(ends[0], ends[1]);
    }

    @NotNull
    @Override
    public Iterator<Pair<FastqRecord, FastqRecord>> iterator() {
        return this;
    }

    @Override
    public boolean hasNext() {
        return nextRecord != null;
    }

    @Override
    public Pair<FastqRecord, FastqRecord> next() {
        Pair<FastqRecord, FastqRecord> p = nextRecord;

        nextRecord = getNextRecord();

        return p;
    }

    public File[] getFile() { return files; }
}
