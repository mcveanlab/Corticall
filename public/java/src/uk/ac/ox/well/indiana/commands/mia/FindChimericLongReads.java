package uk.ac.ox.well.indiana.commands.mia;

import com.google.common.base.Joiner;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class FindChimericLongReads extends Module {
    @Argument(fullName="kmerTable", shortName="t", doc="Kmer table")
    public File KMER_TABLE;

    @Argument(fullName="longReads", shortName="l", doc="Long reads")
    public ArrayList<FastqReader> LONG_READS;

    @Output
    public PrintStream out;

    private Map<CortexKmer, String> loadDiagnosticKmers() {
        Map<CortexKmer, String> kmers = new HashMap<CortexKmer, String>();

        TableReader tr = new TableReader(KMER_TABLE);
        for (Map<String, String> te : tr) {
            CortexKmer ck = new CortexKmer(te.get("kmer"));
            String chr = te.get("chr");

            kmers.put(ck, chr);
        }

        return kmers;
    }

    @Override
    public void execute() {
        log.info("Loading diagnostic kmers...");
        Map<CortexKmer, String> kmers = loadDiagnosticKmers();
        log.info("  {} kmers", kmers.size());

        log.info("Processing long reads...");
        TableWriter tw = new TableWriter(out);

        int kmerSize = kmers.keySet().iterator().next().length();
        for (FastqReader fqr : LONG_READS) {
            log.info("  {}", fqr.getFile().getName());

            Map<Integer, Integer> chrsTouched = new TreeMap<Integer, Integer>();
            int numReads = 0;

            for (FastqRecord fr : fqr) {
                if (numReads % 10000 == 0) {
                    log.info("    {} reads", numReads);
                }
                numReads++;

                Map<String, Integer> chrCounts = new TreeMap<String, Integer>();
                String seq = fr.getReadString();
                for (int i = 0; i <= seq.length() - kmerSize; i++) {
                    CortexKmer ck = new CortexKmer(seq.substring(i, i + kmerSize));

                    if (kmers.containsKey(ck)) {
                        String chr = kmers.get(ck);
                        if (!chrCounts.containsKey(chr)) {
                            chrCounts.put(chr, 1);
                        } else {
                            chrCounts.put(chr, chrCounts.get(chr) + 1);
                        }
                    }
                }

                Map<String, Integer> chrCountsFilt = new TreeMap<String, Integer>();
                for (String chr : chrCounts.keySet()) {
                    if (chrCounts.get(chr) >= 10) {
                        chrCountsFilt.put(chr, chrCounts.get(chr));
                    }
                }

                int numChrs = chrCountsFilt.size();

                if (!chrsTouched.containsKey(numChrs)) {
                    chrsTouched.put(numChrs, 1);
                } else {
                    chrsTouched.put(numChrs, chrsTouched.get(numChrs) + 1);
                }

                //log.info("  {} {} {}", fr.getReadHeader(), seq.length(), Joiner.on(", ").withKeyValueSeparator("=").join(chrCountsFilt));
            }

            log.info("    {}", Joiner.on(", ").withKeyValueSeparator("=").join(chrsTouched));

            Map<String, String> te = new LinkedHashMap<String, String>();
            te.put("filename", fqr.getFile().getName());

            for (int numChrs : chrsTouched.keySet()) {
                te.put(String.valueOf(numChrs), String.valueOf(chrsTouched.get(numChrs)));
            }

            tw.addEntry(te);
        }
    }
}
