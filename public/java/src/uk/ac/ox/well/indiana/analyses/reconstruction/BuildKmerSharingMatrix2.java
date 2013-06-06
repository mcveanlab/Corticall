package uk.ac.ox.well.indiana.analyses.reconstruction;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.io.*;
import java.util.*;

public class BuildKmerSharingMatrix2 extends Tool {
    @Argument(fullName="contigsTable", shortName="ct", doc="Contigs table")
    public File CONTIG_TABLE;

    @Argument(fullName="ignoreNonPanelKmers", shortName="inpk", doc="Only examine sharing for panel kmers.  Ignore non-panel kmers.")
    public Boolean IGNORE_NON_PANEL_KMERS = false;

    @Argument(fullName="samples", shortName="s", doc="Samples")
    public ArrayList<String> SAMPLES;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, Integer> sampleIndex = new TreeMap<String, Integer>();
        Map<Integer, String> indexSample = new HashMap<Integer, String>();
        for (int i = 0; i < SAMPLES.size(); i++) {
            sampleIndex.put(SAMPLES.get(i), i);
            indexSample.put(i, SAMPLES.get(i));
        }

        Map<CortexKmer, boolean[]> kmerPresence = new HashMap<CortexKmer, boolean[]>();

        /*
        try {
            BufferedReader br = new BufferedReader(new FileReader(CONTIG_TABLE));

            String[] header = br.readLine().split("\t");
            String line;
            int numLines = 0;
            while ((line = br.readLine()) != null) {
                if (numLines % 1000 == 0) {
                    log.info("Lines read: {}", numLines);
                }
                numLines++;

                String[] fields = line.split("\t");

                Map<String, String> entry = new HashMap<String, String>();
                for (int i = 0; i < header.length; i++) {
                    entry.put(header[i], fields[i]);
                }

                CortexKmer contig = new CortexKmer(entry.get("contig"));

                String[] kmers = entry.get("kmers").split(",");
                int kmerSize = kmers[0].length();

                String sample = entry.get("sample");
                int sindex = sampleIndex.get(sample);

                for (int i = 0; i < contig.length() - kmerSize; i++) {
                    CortexKmer kmer = contig.getSubKmer(i, kmerSize);

                    if (!kmerPresence.containsKey(kmer)) {
                        kmerPresence.put(kmer, new boolean[SAMPLES.size()]);

                        for (int j = 0; j < SAMPLES.size(); j++) {
                            kmerPresence.get(kmer)[j] = false;
                        }
                    }

                    kmerPresence.get(kmer)[sindex] = true;
                }
            }
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Could not open file '" + CONTIG_TABLE.getAbsolutePath() + "': " + e);
        } catch (IOException e) {
            throw new RuntimeException("Could not read line from '" + CONTIG_TABLE.getAbsolutePath() + "': " + e);
        }
        */

        TableReader tr = new TableReader(CONTIG_TABLE);
        for (Map<String, String> entry : tr) {
            CortexKmer contig = new CortexKmer(entry.get("contig"));

            String[] kmers = entry.get("kmers").split(",");
            int kmerSize = kmers[0].length();

            String sample = entry.get("sample");
            int sindex = sampleIndex.get(sample);

            for (int i = 0; i < contig.length() - kmerSize; i++) {
                CortexKmer kmer = contig.getSubKmer(i, kmerSize);

                if (!kmerPresence.containsKey(kmer)) {
                    kmerPresence.put(kmer, new boolean[SAMPLES.size()]);

                    for (int j = 0; j < SAMPLES.size(); j++) {
                        kmerPresence.get(kmer)[j] = false;
                    }
                }

                kmerPresence.get(kmer)[sindex] = true;
            }
        }

        Map<String, Map<String, Integer>> data = new HashMap<String, Map<String, Integer>>();
        for (String sample1 : SAMPLES) {
            data.put(sample1, new HashMap<String, Integer>());

            for (String sample2 : SAMPLES) {
                data.get(sample1).put(sample2, 0);
            }
        }

        for (CortexKmer kmer : kmerPresence.keySet()) {
            boolean[] presence = kmerPresence.get(kmer);

            for (int i = 0; i < presence.length; i++) {
                String sample1 = indexSample.get(i);

                if (presence[i]) {
                    for (int j = 0; j < presence.length; j++) {
                        String sample2 = indexSample.get(j);

                        if (presence[j]) {
                            int count = data.get(sample1).get(sample2);
                            data.get(sample1).put(sample2, count + 1);
                        }
                    }
                }
            }
        }

        List<String> mheader = new ArrayList<String>();
        mheader.add("\t");
        for (String sample : data.keySet()) {
            mheader.add(sample);
        }

        out.println(Joiner.on("\t").join(mheader));

        for (String sample1 : data.keySet()) {
            List<String> row = new ArrayList<String>();
            row.add(sample1);

            for (String sample2 : data.keySet()) {
                row.add(String.valueOf(data.get(sample1).get(sample2)));
            }

            out.println(Joiner.on("\t").join(row));
        }
    }
}
