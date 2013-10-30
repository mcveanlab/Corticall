package uk.ac.ox.well.indiana.mod;

import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import uk.ac.ox.well.indiana.tools.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.ApproximateHashSet;
import uk.ac.ox.well.indiana.utils.file.FileAndPathUtils;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.performance.PerformanceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class CountKmersApproximately extends Module {
    @Argument(fullName="in", shortName="i", doc="Fasta, fastq, fastq.gz, Cortex binary, or BAM file in which kmers should be counted")
    public ArrayList<File> IN;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Argument(fullName="exact", shortName="e", doc="Enable exact counting")
    public Boolean EXACT = false;

    @Output
    public PrintStream out;

    private class IteratingFileWrapper implements Iterable<String>, Iterator<String> {
        private File file;
        private FileAndPathUtils.FileType ft;
        private String nextRecord;

        private FastaSequenceFile fsr;
        private FastqReader fqr;
        private SAMFileReader sfr;
        private SAMRecordIterator sfrit;
        private CortexGraph cgr;

        public IteratingFileWrapper(File file) {
            this.file = file;

            ft = FileAndPathUtils.getFileType(file);
        }

        @Override
        public Iterator<String> iterator() {
            switch (ft) {
                case FASTA:
                    fsr = new FastaSequenceFile(file, true);
                    break;
                case FASTQ:
                    fqr = new FastqReader(file);
                    break;
                case BAM:
                    sfr = new SAMFileReader(file);
                    sfr.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
                    sfrit = sfr.iterator();
                    break;
                case CTX:
                    cgr = new CortexGraph(file);
                    break;
                default:
                    throw new RuntimeException("Reader for type '" + ft + "' not implemented yet.");
            }

            nextRecord = getNextRecord();

            return this;
        }

        private String getNextRecord() {
            switch (ft) {
                case FASTA:
                    ReferenceSequence seq = fsr.nextSequence();
                    if (seq != null) {
                        return new String(seq.getBases());
                    }
                    break;
                case FASTQ:
                    if (fqr.hasNext()) {
                        FastqRecord fq = fqr.next();
                        return fq.getReadString();
                    }
                    break;
                case BAM:
                    if (sfrit.hasNext()) {
                        SAMRecord sr = sfrit.next();
                        return new String(sr.getReadBases());
                    }
                    break;
                case CTX:
                    if (cgr.hasNext()) {
                        CortexRecord cr = cgr.next();
                        return cr.getKmerAsString();
                    }
                    break;
                default:
                    throw new RuntimeException("Reader for type '" + ft + "' not implemented yet.");
            }

            return null;
        }

        @Override
        public boolean hasNext() {
            return nextRecord != null;
        }

        @Override
        public String next() {
            String currentRecord = nextRecord;

            nextRecord = getNextRecord();
            if (nextRecord == null) {
                close();
            }

            return currentRecord;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }

        public void close() {
            switch (ft) {
                case FASTA:
                    fsr.close();
                    break;
                case FASTQ:
                    fqr.close();
                    break;
                case BAM:
                    sfr.close();
                    break;
                case CTX:
                    cgr.close();
                    break;
                default:
                    throw new RuntimeException("Reader for type '" + ft + "' not implemented yet.");
            }
        }
    }

    @Override
    public void execute() {
        Set<CortexKmer> kmers;
        if (EXACT) {
            kmers = new HashSet<CortexKmer>();
        } else {
            kmers = new ApproximateHashSet<CortexKmer>();
        }

        for (File file : IN) {
            IteratingFileWrapper ifw = new IteratingFileWrapper(file);

            int recordsSeen = 0;
            for (String record : ifw) {
                for (int i = 0; i <= record.length() - KMER_SIZE; i++) {
                    CortexKmer kmer = new CortexKmer(record.substring(i, i + KMER_SIZE));

                    kmers.add(kmer);
                }

                recordsSeen++;

                if (recordsSeen % 1000000 == 0) {
                    log.info("file: {} records: {}, mem: {}", file.getName(), recordsSeen, PerformanceUtils.getCompactMemoryUsageStats());
                }
            }
        }

        out.println("Kmers: " + kmers.size());
    }
}
