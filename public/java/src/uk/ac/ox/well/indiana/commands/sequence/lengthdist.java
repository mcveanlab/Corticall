package uk.ac.ox.well.indiana.commands.sequence;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class lengthdist extends Module {
    @Argument(fullName="file", shortName="f", doc="File (Fastq, or BAM)")
    public HashMap<String, File> FILES;

    @Output
    public PrintStream out;

    @Output
    public PrintStream sout;

    private class GeneralizedSequenceReader implements Iterator<String>, Iterable<String> {
        private SAMFileReader sfr;
        private FastqReader fqr;

        private SAMRecordIterator sri;
        private Iterator<FastqRecord> fri;

        public GeneralizedSequenceReader(File file) {
            if (file.getName().endsWith(".bam")) {
                sfr = new SAMFileReader(file);
                sfr.setValidationStringency(ValidationStringency.SILENT);
            } else if (file.getName().endsWith(".fastq")) {
                fqr = new FastqReader(file, true);
            } else {
                throw new IndianaException("Cannot parse '" + file.getAbsolutePath() + "' with generalized reader");
            }
        }

        @Override
        public Iterator<String> iterator() {
            if (sfr != null) {
                sri = sfr.iterator();
            } else if (fqr != null) {
                fri = fqr.iterator();
            }

            return this;
        }

        @Override
        public boolean hasNext() {
            return (sri != null) ? sri.hasNext() : fri.hasNext();
        }

        @Override
        public String next() {
            return (sri != null) ? sri.next().getReadString() : fri.next().getReadString();
        }

        @Override
        public void remove() {
            throw new IndianaException("Unsupported operation");
        }
    }

    @Override
    public void execute() {
        Map<String, Map<Integer, Integer>> lengthDists = new HashMap<String, Map<Integer, Integer>>();

        TableWriter statsw = new TableWriter(sout);

        int overallMinLength = 0;
        int overallMaxLength = 0;

        log.info("Processing files...");
        for (String key : FILES.keySet()) {
            lengthDists.put(key, new TreeMap<Integer, Integer>());

            GeneralizedSequenceReader gsr = new GeneralizedSequenceReader(FILES.get(key));

            int numReads = 0;
            int minLength = Integer.MAX_VALUE;
            int maxLength = 0;

            Collection<String> reads = new ArrayList<String>();

            for (String read : gsr) {
                int length = read.length();

                if (!lengthDists.get(key).containsKey(length)) {
                    lengthDists.get(key).put(length, 1);
                } else {
                    lengthDists.get(key).put(length, lengthDists.get(key).get(length) + 1);
                }

                if (length < minLength) { minLength = length; }
                if (length > maxLength) { maxLength = length; }

                if (length < overallMinLength) { overallMinLength = length; }
                if (length > overallMaxLength) { overallMaxLength = length; }

                numReads++;

                reads.add(read);
            }

            float meanLength = SequenceUtils.meanLength(reads);
            int n50Length = SequenceUtils.computeN50Length(reads);
            int n50Value = SequenceUtils.computeN50Value(reads);

            Map<String, String> te = new LinkedHashMap<String, String>();
            te.put("key", key);
            te.put("numReads", String.valueOf(numReads));
            te.put("minLength", String.valueOf(minLength));
            te.put("maxLength", String.valueOf(maxLength));
            te.put("meanLength", String.valueOf(meanLength));
            te.put("n50Length", String.valueOf(n50Length));
            te.put("n50Value", String.valueOf(n50Value));

            statsw.addEntry(te);
        }

        TableWriter histw = new TableWriter(out);

        for (int i = overallMinLength; i <= overallMaxLength; i++) {
            Map<String, String> te = new LinkedHashMap<String, String>();

            te.put("length", String.valueOf(i));

            for (String key : lengthDists.keySet()) {
                te.put(key, lengthDists.get(key).containsKey(i) ? String.valueOf(lengthDists.get(key).get(i)) : "0");
            }

            histw.addEntry(te);
        }
    }
}
