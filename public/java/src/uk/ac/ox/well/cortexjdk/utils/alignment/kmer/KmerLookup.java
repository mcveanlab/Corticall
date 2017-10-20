package uk.ac.ox.well.cortexjdk.utils.alignment.kmer;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.cortexjdk.utils.alignment.pairwise.BwaAligner;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.utils.LineReader;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class KmerLookup {
    private File refFile;
    private IndexedFastaSequenceFile ref;
    private Set<String> sources = new HashSet<>();
    private BwaAligner bwa;

    public KmerLookup(File refFile) { initialize(refFile); }

    private void initialize(File refFile) {
        try {
            this.refFile = refFile;
            this.ref = new IndexedFastaSequenceFile(refFile);

            if (this.ref.getSequenceDictionary() == null) {
                throw new CortexJDKException("Reference must have a sequence dictionary");
            }

            File sourcesFile = new File(refFile.getAbsolutePath() + ".sources");
            if (!sourcesFile.exists()) {
                throw new CortexJDKException("Fasta sequence must be indexed with IndexReference first");
            } else {
                LineReader lr = new LineReader(sourcesFile);
                while (lr.hasNext()) {
                    String line = lr.getNextRecord();

                    sources.add(line);
                }
            }

            bwa = new BwaAligner(refFile.getAbsolutePath());
        } catch (FileNotFoundException e) {
            throw new CortexJDKException("File not found", e);
        }
    }

    public Set<String> getSources() { return sources; }

    public IndexedFastaSequenceFile getReferenceSequence() {
        return ref;
    }

    public BwaAligner getAligner() { return bwa; }

    public String find(Interval interval) {
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

    public Set<Interval> find(String seq) {
        List<SAMRecord> alignments = bwa.align(seq);
        Set<Interval> intervals = new HashSet<>();

        for (SAMRecord sr : alignments) {
            if (sr.getIntegerAttribute("NM") == 0 && sr.getCigar().numCigarElements() == 1) {
                intervals.add(new Interval(sr.getContig(), sr.getAlignmentStart() + 1, sr.getAlignmentEnd() + 1, sr.getReadNegativeStrandFlag(), null));
            }
        }

        return intervals;
    }

    public static File createIndex(File refFile, String... sources) {
        try {
            File sourcesFile = new File(refFile.getAbsolutePath() + ".sources");
            PrintStream ps = new PrintStream(sourcesFile);

            for (String source : sources) {
                ps.println(source);
            }

            ps.close();

            return sourcesFile;
        } catch (FileNotFoundException e) {
            throw new CortexJDKException("File not found", e);
        }
    }
}
