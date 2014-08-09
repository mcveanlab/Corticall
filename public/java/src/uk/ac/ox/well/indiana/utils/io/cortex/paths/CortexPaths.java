package uk.ac.ox.well.indiana.utils.io.cortex.paths;

import org.json.JSONArray;
import org.json.JSONObject;
import uk.ac.ox.well.indiana.Indiana;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexColor;

import java.io.*;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPInputStream;

public class CortexPaths implements Iterable<CortexPathsRecord>, Iterator<CortexPathsRecord> {
    private File pathsFile;

    private int version;
    private int kmerSize;
    private int numColors;
    private long numKmersInGraph;
    private long numKmersWithPaths;
    private long numPaths;
    private long pathBytes;

    private List<CortexColor> colors = new ArrayList<CortexColor>();

    private BufferedReader buffered;

    private CortexPathsRecord nextRecord = null;
    private int recordsSeen = 0;
    private long recordsStart = 0;

    public CortexPaths(String pathsString) {
        this.pathsFile = new File(pathsString);
        loadCortexPaths(this.pathsFile);
    }

    public CortexPaths(File pathsFile) {
        this.pathsFile = pathsFile;
        loadCortexPaths(this.pathsFile);
    }

    private void loadCortexPaths(File cortexPaths) {
        try {
            InputStream fileStream = new FileInputStream(cortexPaths);
            InputStream gzipStream = new GZIPInputStream(fileStream);
            Reader decoder = new InputStreamReader(gzipStream, Charset.forName("UTF-8"));
            buffered = new BufferedReader(decoder);

            StringBuilder headerBuilder = new StringBuilder();
            boolean inHeader = false;

            String line;
            while ((line = buffered.readLine()) != null) {
                recordsStart += line.length() + 1;

                if (line.equals("{")) { inHeader = true; }
                if (inHeader) { headerBuilder.append(line); }
                if (line.equals("}")) { break; }
            }

            JSONObject header = new JSONObject(headerBuilder.toString());

            version = header.getInt("formatVersion");

            if (version != 2) {
                throw new IndianaException("Cannot parse CortexPaths format version '" + version + "'");
            }

            numColors = header.getInt("ncols");
            kmerSize = header.getInt("kmer_size");
            numKmersInGraph = header.getLong("num_kmers_in_graph");
            numKmersWithPaths = header.getLong("num_kmers_with_paths");
            numPaths = header.getLong("num_paths");
            pathBytes = header.getLong("path_bytes");

            JSONArray colors = header.getJSONArray("colours");
            for (int i = 0; i < colors.length(); i++) {
                JSONObject jColor = colors.getJSONObject(i);

                CortexColor color = new CortexColor();
                color.setSampleName(jColor.getString("sample"));
                color.setTotalSequence(jColor.getLong("total_sequence"));
                color.setTopClippingApplied(true);
                color.setLowCovgSupernodesRemoved(true);

                this.colors.add(color);
            }

            while ((line = buffered.readLine()) != null) {
                if (!line.isEmpty() && !line.startsWith("#")) {
                    break;
                }

                recordsStart += line.length() + 1;
                buffered.mark(100);
            }
        } catch (IOException e) {
            throw new IndianaException("Unable to load Cortex paths file '" + pathsFile.getAbsolutePath() + "'", e);
        }
    }

    private void moveToBeginningOfRecordsSection() {
        try {
            buffered.reset();
            recordsSeen = 0;
        } catch (IOException e) {
            try {
                buffered.close();

                InputStream fileStream = new FileInputStream(this.pathsFile);
                InputStream gzipStream = new GZIPInputStream(fileStream);
                Reader decoder = new InputStreamReader(gzipStream, Charset.forName("UTF-8"));
                buffered = new BufferedReader(decoder);

                buffered.skip(recordsStart);
                recordsSeen = 0;
            } catch (IOException e2) {
                throw new IndianaException("Unable to restart iteration over CortexPaths file", e2);
            }
        }

        nextRecord = getNextRecord();
    }

    private CortexPathsRecord getNextRecord() {
        if (recordsSeen < numKmersWithPaths) {
            try {
                String[] kmerLine = buffered.readLine().split("\\s+");

                String kmer = kmerLine[0];
                int numPaths = Integer.valueOf(kmerLine[1]);

                List<CortexJunctionsRecord> cjs = new ArrayList<CortexJunctionsRecord>();

                for (int i = 0; i < numPaths; i++) {
                    String[] pathLine = buffered.readLine().split("\\s+");

                    String orientation = pathLine[0];
                    int numKmers = Integer.valueOf(pathLine[1]);
                    int numJunctions = Integer.valueOf(pathLine[2]);
                    int[] coverages = new int[numColors];

                    for (int c = 0; c < numColors; c++) {
                        coverages[c] = Integer.valueOf(pathLine[3 + c]);
                    }

                    String junctions = pathLine[pathLine.length - 1];

                    CortexJunctionsRecord cj = new CortexJunctionsRecord(orientation.equals("F"), numKmers, numJunctions, coverages, junctions);
                    cjs.add(cj);
                }

                recordsSeen++;

                return new CortexPathsRecord(kmer, cjs);
            } catch (IOException e) {
                throw new IndianaException("Unable to parse CortexPaths record", e);
            }
        }

        return null;
    }

    @Override
    public Iterator<CortexPathsRecord> iterator() {
        moveToBeginningOfRecordsSection();

        return this;
    }

    @Override
    public boolean hasNext() {
        return nextRecord != null;
    }

    @Override
    public CortexPathsRecord next() {
        CortexPathsRecord currentRecord = nextRecord;

        nextRecord = getNextRecord();
        if (nextRecord == null) {
            close();
        }

        return currentRecord;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    public void close() {
        try {
            this.buffered.close();
        } catch (IOException e) {
            throw new IndianaException("Unable to close CortexPaths file", e);
        }
    }

    public String toString() {
        StringBuilder out = new StringBuilder();

        out.append("file: ").append(pathsFile.getAbsolutePath()).append("\n");
        out.append("----\n");
        out.append("binary version: ").append(version).append("\n");
        out.append("kmer size: ").append(kmerSize).append("\n");
        out.append("path bytes: ").append(pathBytes).append("\n");
        out.append("colors: ").append(numColors).append("\n");
        for (int i = 0; i < colors.size(); i++) {
            CortexColor color = colors.get(i);

            out.append("-- Color ").append(i).append(" --\n");
            out.append("  sample: '").append(color.getSampleName()).append("'\n");
            out.append("  total_sequence: ").append(color.getTotalSequence()).append("\n");
            out.append("  cleaned_tips: ").append(color.isTopClippingApplied()).append("\n");
            out.append("  cleaned_supernodes: ").append(color.isLowCovgSupernodesRemoved()).append("\n");
        }
        out.append("----\n");
        out.append("kmers in graph: ").append(numKmersInGraph).append("\n");
        out.append("kmers with paths: ").append(numKmersWithPaths).append("\n");
        out.append("paths: ").append(numPaths).append("\n");
        out.append("----\n");

        return out.toString();
    }

    public int getVersion() { return version; }
    public int getKmerSize() { return kmerSize; }
    public int getNumColors() { return numColors; }
    public long getNumKmersInGraph() { return numKmersInGraph; }
    public long getNumKmersWithPaths() { return numKmersWithPaths; }
    public long getNumPaths() { return numPaths; }
    public long getPathBytes() { return pathBytes; }
    public File getCortexPathsFile() { return pathsFile; }

    public boolean hasColor(int color) {
        return (color < colors.size());
    }

    public CortexColor getColor(int color) {
        return colors.get(color);
    }

    public int getColorForSampleName(String sampleName) {
        int sampleColor = -1;
        int sampleCopies = 0;

        for (int color = 0; color < colors.size(); color++) {
            if (colors.get(color).getSampleName().equalsIgnoreCase(sampleName)) {
                sampleColor = color;
                sampleCopies++;
            }
        }

        return (sampleCopies == 1) ? sampleColor : -1;
    }
}
