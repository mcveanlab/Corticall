package uk.ac.ox.well.indiana.utils.io.cortex;

import uk.ac.ox.well.indiana.utils.io.utils.BinaryFile;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

public class CortexGraph { //implements Iterable<CortexRecord> {
    private File cortexFile;
    private BinaryFile in;

    private int version;
    private int kmerSize;
    private int kmerBits;
    private int numColors;
    private ArrayList<CortexColor> colors = new ArrayList<CortexColor>();
    private long dataOffset;

    public CortexGraph(String cortexFilePath) {
        this.cortexFile = new File(cortexFilePath);
        loadCortexGraph(this.cortexFile);
    }

    public CortexGraph(File cortexFile) {
        this.cortexFile = cortexFile;
        loadCortexGraph(this.cortexFile);
    }

    private void loadCortexGraph(File cortexFile) {
        try {
            in = new BinaryFile(cortexFile, "r");

            byte[] headerStart = new byte[6];
            in.read(headerStart);
            String headerStartStr = new String(headerStart);

            if (!headerStartStr.equalsIgnoreCase("CORTEX")) {
                throw new RuntimeException("The file '" + cortexFile.getAbsolutePath() + "' does not appear to be a Cortex graph");
            }

            version = in.readUnsignedInt();
            kmerSize = in.readUnsignedInt();
            kmerBits = in.readUnsignedInt();
            numColors = in.readUnsignedInt();

            for (int color = 0; color < numColors; color++) {
                colors.add(new CortexColor());
            }

            for (int color = 0; color < numColors; color++) {
                colors.get(color).setMeanReadLength(in.readUnsignedInt());
            }

            for (int color = 0; color < numColors; color++) {
                colors.get(color).setTotalSequence(in.readUnsignedLong());
            }

            for (int color = 0; color < numColors; color++) {
                int sampleNameLength = in.readUnsignedInt();
                byte[] sampleName = new byte[sampleNameLength];
                in.read(sampleName);

                colors.get(color).setSampleName(new String(sampleName));
            }

            // Todo: fix this at some point - we're not actually getting the error rate properly
            for (int color = 0; color < numColors; color++) {
                byte[] errorRate = new byte[16];
                in.read(errorRate);
            }

            for (int color = 0; color < numColors; color++) {
                colors.get(color).setTopClippingApplied(in.readBoolean());
                colors.get(color).setLowCovgSupernodesRemoved(in.readBoolean());
                colors.get(color).setLowCovgKmersRemoved(in.readBoolean());
                colors.get(color).setCleanedAgainstGraph(in.readBoolean());
                colors.get(color).setLowCovSupernodesThreshold(in.readUnsignedInt());
                colors.get(color).setLowCovKmerThreshold(in.readUnsignedInt());

                int graphNameLength = in.readUnsignedInt();
                byte[] graphName = new byte[graphNameLength];
                in.read(graphName);

                colors.get(color).setCleanedAgainstGraphName(new String(graphName));
            }

            byte[] headerEnd = new byte[6];
            in.read(headerEnd);
            String headerEndStr = new String(headerEnd);

            if (!headerEndStr.equalsIgnoreCase("CORTEX")) {
                throw new RuntimeException("We didn't see a proper header terminator at the expected place in Cortex graph '" + cortexFile.getAbsolutePath() + "'");
            }

            dataOffset = in.getFilePointer();
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Cortex graph file '" + cortexFile.getAbsolutePath() + "' not found: " + e);
        } catch (IOException e) {
            throw new RuntimeException("Error while parsing Cortex graph file '" + cortexFile.getAbsolutePath() + "': " + e);
        }
    }

    public File getCortexFile() {
        return cortexFile;
    }

    public int getVersion() {
        return version;
    }

    public int getKmerSize() {
        return kmerSize;
    }

    public int getKmerBits() {
        return kmerBits;
    }

    public int getNumColors() {
        return numColors;
    }

    public ArrayList<CortexColor> getColors() {
        return colors;
    }

    public long getDataOffset() {
        return dataOffset;
    }

    public String toString() {
        String info = "file: " + cortexFile.getAbsolutePath() + "\n"
                    + "----" + "\n"
                    + "binary version: " + this.getVersion() + "\n"
                    + "kmer size: " + this.getKmerSize() + "\n"
                    + "bitfields: " + this.getKmerBits() + "\n"
                    + "colors: " + this.getNumColors() + "\n";

        for (int color = 0; color < this.getNumColors(); color++) {
            CortexColor cortexColor = this.getColors().get(color);
            info += "-- Color " + color + " --\n"
                 +  "  sample name: '" + cortexColor.getSampleName() + "'\n"
                 +  "  mean read length: " + cortexColor.getMeanReadLength() + "\n"
                 +  "  total sequence loaded: " + cortexColor.getTotalSequence() + "\n"
                 +  "  sequence error rate: " + "(not parsed)" + "\n"
                 +  "  tip clipping: " + (cortexColor.isTopClippingApplied() ? "yes" : "no") + "\n"
                 +  "  remove_low_coverage_supernodes: " + (cortexColor.isLowCovgSupernodesRemoved() ? "yes" : "no") + "\n"
                 +  "  remove_low_coverage_kmers: " + (cortexColor.isLowCovgKmersRemoved() ? "yes" : "no") + "\n"
                 +  "  cleaned against graph: " + (cortexColor.isCleanedAgainstGraph() ? "yes" : "no") + "\n";
        }

        info += "----" + "\n";

        return info;
    }

    public CortexRecord next() {
        try {
            long[] binaryKmer = new long[kmerBits];
            for (int bits = 0; bits < kmerBits; bits++) {
                binaryKmer[bits] = in.readUnsignedLong();
            }

            int[] coverages = new int[numColors];
            for (int color = 0; color < numColors; color++) {
                coverages[color] = in.readUnsignedInt();
            }

            byte[] edges = new byte[numColors];
            for (int color = 0; color < numColors; color++) {
                edges[color] = (byte) in.readUnsignedByte();
            }

            return new CortexRecord(binaryKmer, coverages, edges, kmerSize, kmerBits);

            /*
            byte[] kmer = new byte[kmerSize];
            for (int i = kmerSize - 1; i >= 0; i--) {
                kmer[i] = binaryNucleotideToChar(binaryKmer[kmerBits - 1] & 0x3);
                shiftBinaryKmerByOneBase(binaryKmer, kmerBits);
            }
            String kmerString = new String(kmer);

            byte[] str = {'a', 'c', 'g', 't', 'A', 'C', 'G', 'T'};
            String[] edges = new String[numColors];
            for (int color = 0; color < numColors; color++) {
                byte edge = (byte) in.readUnsignedByte();

                int left = (edge >> 4);
                int right = (edge & 0xf);

                byte[] edgeStr = new byte[8];

                for (int i = 0; i < 4; i++) {
                    int leftEdge = (left & (0x1 << (3-i)));
                    edgeStr[i] = (byte) ((leftEdge != 0) ? str[i] : '.');

                    int rightEdge = (right & (0x1 << i));
                    edgeStr[i+4] = (byte) ((rightEdge != 0) ? str[i+4] : '.');
                }

                edges[color] = new String(edgeStr);
            }

            System.out.println("Record:");
            System.out.println("\t" + kmerString);
            for (int color = 0; color < numColors; color++) {
                System.out.println("\t" + coverages[color]);
                System.out.println("\t" + edges[color]);
            }
            */
        } catch (IOException e) {
            throw new RuntimeException("Error while parsing Cortex record");
        }
    }
}
