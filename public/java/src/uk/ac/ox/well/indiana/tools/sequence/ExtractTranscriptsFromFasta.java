package uk.ac.ox.well.indiana.tools.sequence;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class ExtractTranscriptsFromFasta extends Tool {
    @Argument(fullName="fasta", shortName="f", doc="Fasta file from which sequences should be extracted")
    public FastaSequenceFile FASTA;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE;

    @Argument(fullName="gff", shortName="gff", doc="The GFF file with the transcript definitions")
    public File REGIONS;

    @Output
    public PrintStream out;

    private class Transcript {
        String name;
        String contig;
        int start;
        int stop;

        public Transcript(String name, String contig, int start, int stop) {
            this.name = name;
            this.contig = contig;
            this.start = start;
            this.stop = stop;
        }

        public String getName() { return name; }

        public String getContig() { return contig; }

        public int getStart() { return start; }

        public int getStop() { return stop; }

        public String toString() {
            return getName() + " " + getContig() + ":" + getStart() + "-" + getStop();
        }
    }

    private HashMap<String, ArrayList<Transcript>> loadGFFFile(File regions) {
        HashMap<String, ArrayList<Transcript>> transcripts = new HashMap<String, ArrayList<Transcript>>();

        try {
            BufferedReader gffReader = new BufferedReader(new FileReader(regions));

            String line;
            while ((line = gffReader.readLine()) != null) {
                if (!line.startsWith("#") && line.contains("gene")) {
                    String[] fields = line.split("\\s+");
                    String contigName = fields[0];
                    int start = Integer.valueOf(fields[3]);
                    int stop = Integer.valueOf(fields[4]);

                    String transcriptName = "unknown";

                    String[] info = fields[8].split(";");
                    for (String infoField : info) {
                        if (infoField.startsWith("ID=")) {
                            String[] infoFieldPieces = infoField.split("=");

                            transcriptName = infoFieldPieces[1];
                        }
                    }

                    Transcript t = new Transcript(transcriptName, contigName, start, stop);

                    if (!transcripts.containsKey(contigName)) {
                        transcripts.put(contigName, new ArrayList<Transcript>());
                    }

                    transcripts.get(contigName).add(t);
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return transcripts;
    }

    @Override
    public int execute() {
        HashMap<String, ArrayList<Transcript>> transcripts = loadGFFFile(REGIONS);

        ArrayList<ReferenceSequence> contigs = new ArrayList<ReferenceSequence>();
        HashMap<Integer, Integer> kmerHash = new HashMap<Integer, Integer>();

        ReferenceSequence c;
        while ((c = FASTA.nextSequence()) != null) {
            System.out.println(c);
            contigs.add(c);

            byte[] b = c.getBases();
            for (int i = 0; i < b.length - KMER_SIZE; i++) {
                byte[] b1 = Arrays.copyOfRange(b, i, i + KMER_SIZE);

                String b2 = new String(SequenceUtils.getAlphanumericallyLowestOrientation(b1));
                int b2hashCode = b2.hashCode();

                if (!kmerHash.containsKey(b2hashCode)) {
                    kmerHash.put(b2hashCode, 1);
                } else {
                    kmerHash.put(b2hashCode, kmerHash.get(b2hashCode) + 1);
                }
            }
        }

        for (ReferenceSequence contig : contigs) {
            System.out.println(contig);

            String[] name = contig.getName().split("\\s+");
            String contigName = name[0];

            byte[] bases = contig.getBases();

            if (transcripts.containsKey(contigName)) {
                for (Transcript t : transcripts.get(contigName)) {
                    byte[] transcriptBases = Arrays.copyOfRange(bases, t.getStart()-1, t.getStop());

                    for (int i = 0; i < transcriptBases.length - KMER_SIZE; i++) {
                        String kmer = new String(SequenceUtils.getAlphanumericallyLowestOrientation(Arrays.copyOfRange(transcriptBases, i, i + KMER_SIZE)));

                        if (kmerHash.containsKey(kmer.hashCode()) && kmerHash.get(kmer.hashCode()) == 1) {
                            out.println(">" + t.getName() + "." + kmer.hashCode());
                            out.println(kmer);
                        }
                    }
                }
            }
        }

        return 0;
    }
}
