package uk.ac.ox.well.indiana.attic.analyses.kmerSharing;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.io.File;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class FindSupernodesWithPanelKmers extends Module {
    @Argument(fullName="supernodes", shortName="s", doc="Supernodes FASTA file")
    public FastaSequenceFile SUPERNODES;

    @Argument(fullName="kmerReferencePanel", shortName="krp", doc="Kmer reference panel")
    public File KMER_REFERENCE_PANEL;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableReader tr = new TableReader(KMER_REFERENCE_PANEL);
        Set<CortexKmer> kmers = new HashSet<CortexKmer>();

        int kmerSize = 0;
        for (Map<String, String> te : tr) {
            CortexKmer kmer = new CortexKmer(te.get("kmer"));

            kmers.add(kmer);

            if (kmerSize == 0) {
                kmerSize = kmer.length();
            }
        }

        ReferenceSequence rseq;
        while ((rseq = SUPERNODES.nextSequence()) != null) {
            String seq = new String(rseq.getBases());

            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + kmerSize));

                if (kmers.contains(kmer)) {
                    out.println(">" + rseq.getName() + "\n" + seq);
                    break;
                }
            }
        }
    }
}
