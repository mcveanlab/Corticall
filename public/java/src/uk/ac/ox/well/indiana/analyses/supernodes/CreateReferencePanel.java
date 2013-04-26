package uk.ac.ox.well.indiana.analyses.supernodes;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

public class CreateReferencePanel extends Tool {
    @Argument(fullName="reference", shortName="R", doc="Reference genome")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="gff", shortName="gff", doc="GFF")
    public GFF3 GFF;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 31;

    @Argument(fullName="padding", shortName="p", doc="Amount to pad the sequences")
    public Integer PADDING = 0;

    @Output
    public PrintStream out;

    private Map<String, Integer> loadKmers() {
        Map<String, Integer> kmers = new HashMap<String, Integer>();

        ReferenceSequence ref;
        while ((ref = REFERENCE.nextSequence()) != null) {
            log.info("Loading kmers for '{}'", ref.getName());

            String seq = new String(ref.getBases());

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                String fw = SequenceUtils.alphanumericallyLowestOrientation(seq.substring(i, i + KMER_SIZE));

                Integer coverage = (kmers.containsKey(fw) ? kmers.get(fw) : 0);
                coverage++;

                kmers.put(fw, coverage);
            }
        }

        return kmers;
    }

    @Override
    public void execute() {
        Map<String, Integer> kmers = loadKmers();

        for (GFF3Record r : GFF) {
            if ("gene".equalsIgnoreCase(r.getType())) {
                String id = r.getAttribute("ID");

                log.info("Processing gene '{}'", id);

                String seq = new String(REFERENCE.getSubsequenceAt(r.getSeqid(), r.getStart(), r.getEnd()).getBases());

                String padding = "";
                for (int i = 0; i < PADDING; i++) {
                    padding += 'N';
                }

                out.println(">" + id);
                out.println(padding + seq + padding);

                for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                    String fw = SequenceUtils.alphanumericallyLowestOrientation(seq.substring(i, i + KMER_SIZE));

                    if (kmers.containsKey(fw) && kmers.get(fw) == 1) {
                        out.println(">" + id + "." + fw);
                        out.println(fw);
                    }
                }
            }
        }
    }
}
