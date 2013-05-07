package uk.ac.ox.well.indiana.analyses.reconstruction;

import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.utils.TableReader2;

import java.io.File;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class BuildKmerSharingMatrix extends Tool {
    @Argument(fullName="contigsTable", shortName="ct", doc="Contigs table")
    public Collection<File> CONTIG_FILES;

    @Argument(fullName="kmerReferencePanel", shortName="krp", doc="Kmer reference panel")
    public File KMER_REFERENCE_PANEL;

    @Argument(fullName="ignoreNonPanelKmers", shortName="inpk", doc="Only examine sharing for panel kmers.  Ignore non-panel kmers.")
    public Boolean USE_NON_PANEL_KMERS = true;

    public Set<CortexKmer> loadKmerReferencePanel() {
        Set<CortexKmer> krp = new HashSet<CortexKmer>();

        TableReader2 tr = new TableReader2(KMER_REFERENCE_PANEL);

        for (Map<String, String> te : tr) {
            CortexKmer ck = new CortexKmer(te.get("kmer"));

            krp.add(ck);
        }

        return krp;
    }

    @Override
    public void execute() {
        Set<CortexKmer> krp = loadKmerReferencePanel();
        int kmerSize = krp.iterator().next().length();

        for (File contigFile : CONTIG_FILES) {
            TableReader2 tr = new TableReader2(contigFile);

            for (Map<String, String> te : tr) {
                CortexKmer ck = new CortexKmer(te.get("contig"));
                int color = Integer.valueOf(te.get("color"));

                for (int i = 0; i <= ck.length() - kmerSize; i++) {
                    CortexKmer sk = ck.getSubKmer(i, kmerSize);

                    if (krp.contains(sk) || USE_NON_PANEL_KMERS) {

                    }
                }
            }
        }
    }
}
