package uk.ac.ox.well.indiana.attic.analyses.reconstruction;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class BuildKmerSharingMatrix extends Module {
    @Argument(fullName="contigsTable", shortName="ct", doc="Contigs table")
    public File CONTIG_TABLE;

    @Argument(fullName="ignoreNonPanelKmers", shortName="inpk", doc="Only examine sharing for panel kmers.  Ignore non-panel kmers.")
    public Boolean IGNORE_NON_PANEL_KMERS = false;

    @Output
    public PrintStream out;

    private class KmerInfo {
        int count = 0;
        boolean isPanelKmer;
        boolean[] presenceInSample;
        String genes;
    }

    private Set<CortexKmer> getCommaDelimitedKmersAsCortexKmers(String commaDelimitedKmers) {
        String[] kmers = commaDelimitedKmers.split(",");

        Set<CortexKmer> ckmers = new HashSet<CortexKmer>();
        for (String kmer : kmers) {
            ckmers.add(new CortexKmer(kmer));
        }

        return ckmers;
    }

    private Set<String> getCommaDelimitedGenes(String commaDelimitedGenes) {
        String[] genes = commaDelimitedGenes.split(",");

        Set<String> cgenes = new TreeSet<String>();
        cgenes.addAll(Arrays.asList(genes));

        return cgenes;
    }

    public void execute() {
        Map<CortexKmer, KmerInfo> kmerInfo = new HashMap<CortexKmer, KmerInfo>();
        Set<String> sampleSet = new TreeSet<String>();
        List<String> samples = new ArrayList<String>();

        TableReader tr = new TableReader(CONTIG_TABLE);

        for (int pass = 0; pass <= 1; pass++) {
            int recIndex = 0;

            samples = (pass == 1) ? new ArrayList<String>(sampleSet) : null;

            for (Map<String, String> te : tr) {
                if (recIndex % (tr.size()/5) == 0) {
                    log.info("Pass {}/2: processed {}/{} records", pass+1, recIndex, tr.size());
                }

                CortexKmer ck = new CortexKmer(te.get("contig"));
                Set<CortexKmer> panelKmers = getCommaDelimitedKmersAsCortexKmers(te.get("kmers"));
                Set<String> genes = getCommaDelimitedGenes(te.get("genes"));
                String sample = te.get("sample");

                int kmerSize = panelKmers.iterator().next().length();

                for (int i = 0; i <= ck.length() - kmerSize; i++) {
                    CortexKmer kmer = ck.getSubKmer(i, kmerSize);

                    if (panelKmers.contains(kmer) || !IGNORE_NON_PANEL_KMERS) {
                        KmerInfo ki = kmerInfo.containsKey(kmer) ? kmerInfo.get(kmer) : new KmerInfo();

                        if (pass == 0) {
                            ki.count++;
                            ki.isPanelKmer = panelKmers.contains(kmer);
                            ki.genes = Joiner.on(",").join(genes);

                            kmerInfo.put(kmer, ki);
                            sampleSet.add(sample);
                        } else {
                            if (sampleSet.contains(sample) && samples != null) {
                                if (ki.presenceInSample == null) {
                                    ki.presenceInSample = new boolean[samples.size()];

                                    for (int index = 0; index < samples.size(); index++) {
                                        ki.presenceInSample[index] = false;
                                    }
                                }

                                int index = samples.indexOf(sample);
                                ki.presenceInSample[index] = true;
                            } else {
                                throw new RuntimeException("Unable to find index for sample '" + sample + "'");
                            }
                        }
                    }
                }

                recIndex++;
            }
        }

        log.info("Writing kmer sharing matrix to disk...");

        TableWriter tw = new TableWriter(out);

        for (CortexKmer ck : kmerInfo.keySet()) {
            KmerInfo ki = kmerInfo.get(ck);

            Map<String, String> te = new LinkedHashMap<String, String>();
            te.put("kmer", ck.getKmerAsString());
            te.put("isPanelKmer", ki.isPanelKmer ? "true" : "false");
            te.put("genes", ki.genes);

            for (int index = 0; index < samples.size(); index++) {
                String sample = samples.get(index);

                te.put(sample, ki.presenceInSample[index] ? "1" : "0");
            }

            tw.addEntry(te);
        }
    }
}
