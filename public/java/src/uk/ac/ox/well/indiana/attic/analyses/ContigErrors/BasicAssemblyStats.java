package uk.ac.ox.well.indiana.attic.analyses.ContigErrors;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

public class BasicAssemblyStats extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (FASTA)")
    public TreeMap<String, FastaSequenceFile> CONTIGS;

    @Argument(fullName="lengthMin", shortName="lmin", doc="Minimum contig length")
    public Integer LENGTH_MIN = 0;

    @Argument(fullName="reference", shortName="r", doc="Reference (FASTA)", required=false)
    public FastaSequenceFile REFERENCE;

    @Output
    public PrintStream out;

    private long getReferenceLength() {
        long totalLength = 0;

        if (REFERENCE != null) {
            ReferenceSequence rseq;
            while ((rseq = REFERENCE.nextSequence()) != null) {
                totalLength += rseq.length();
            }
        }

        return totalLength;
    }

    @Override
    public void execute() {
        TableWriter tw = new TableWriter(out);

        for (String id : CONTIGS.keySet()) {
            Map<String, String> entry = new LinkedHashMap<String, String>();

            List<String> seqs = new ArrayList<String>();

            FastaSequenceFile contigs = CONTIGS.get(id);
            ReferenceSequence rseq;
            long totalSequence = 0;
            while ((rseq = contigs.nextSequence()) != null) {
                String seq = new String(rseq.getBases());

                if (seq.length() >= LENGTH_MIN) {
                    seqs.add(seq);

                    totalSequence += seq.length();
                }
            }

            int numContigs = seqs.size();
            int minLength = SequenceUtils.minLength(seqs);
            int maxLength = SequenceUtils.maxLength(seqs);
            float meanLength = SequenceUtils.meanLength(seqs);
            int n50 = SequenceUtils.computeN50Value(seqs);
            int ng50 = (REFERENCE != null) ? SequenceUtils.computeNG50Value(seqs, getReferenceLength()) : 0;

            entry.put("id", id);
            entry.put("numContigs", String.valueOf(numContigs));
            entry.put("minLength", String.valueOf(minLength));
            entry.put("maxLength", String.valueOf(maxLength));
            entry.put("meanLength", String.format("%.2f", meanLength));
            entry.put("n50", String.valueOf(n50));
            entry.put("ng50", String.valueOf(ng50));
            entry.put("totalSequence", String.valueOf(totalSequence));

            tw.addEntry(entry);
        }
    }
}
