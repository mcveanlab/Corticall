package uk.ac.ox.well.indiana.utils.io.cortex.links;

import org.json.JSONArray;
import org.json.JSONObject;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexColor;

import java.io.*;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPInputStream;

public class CortexLinks implements Iterable<CortexLinksRecord>, Iterator<CortexLinksRecord> {
    private File linksFile;

    private int version;
    private int kmerSize;
    private int numColors;
    private long numKmersInGraph;
    private long numKmersWithLinks;
    private long numLinks;
    private long linkBytes;

    private List<CortexColor> colors = new ArrayList<CortexColor>();

    private BufferedReader buffered;

    private CortexLinksRecord nextRecord = null;
    private int recordsSeen = 0;
    private long recordsStart = 0;

    public CortexLinks(String linksString) {
        this.linksFile = new File(linksString);
        loadCortexLinks(this.linksFile);
    }

    public CortexLinks(File linksFile) {
        this.linksFile = linksFile;
        loadCortexLinks(this.linksFile);
    }

    private void loadCortexLinks(File cortexLinks) {
        try {
            InputStream fileStream = new FileInputStream(cortexLinks);
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

            if (!header.has("formatVersion") && !header.has("format_version")) {
                throw new IndianaException("Cannot parse CortexLinks format version field");
            }
            version = header.has("formatVersion") ? header.getInt("formatVersion") : header.getInt("format_version");

            if (version != 2 && version != 3) {
                throw new IndianaException("Cannot parse CortexLinks format version '" + version + "'");
            }

            if (version == 2) {
                numColors = header.getInt("ncols");
                kmerSize = header.getInt("kmer_size");
                numKmersInGraph = header.getLong("num_kmers_in_graph");
                numKmersWithLinks = header.getLong("num_kmers_with_paths");
                numLinks = header.getLong("num_paths");
                linkBytes = header.getLong("path_bytes");

                JSONArray colors = header.getJSONArray("colours");
                for (int i = 0; i < colors.length(); i++) {
                    JSONObject jColor = colors.getJSONObject(i);

                    CortexColor color = new CortexColor();
                    color.setSampleName(jColor.getString("sample"));
                    color.setTotalSequence(jColor.getLong("total_sequence"));
                    color.setTipClippingApplied(true);
                    color.setLowCovgSupernodesRemoved(true);

                    this.colors.add(color);
                }
            } else {
                JSONObject graph = header.getJSONObject("graph");
                JSONObject paths = header.getJSONObject("paths");

                numColors = graph.getInt("num_colours");
                kmerSize = graph.getInt("kmer_size");
                numKmersInGraph = graph.getLong("num_kmers_in_graph");
                numKmersWithLinks = paths.getLong("num_kmers_with_paths");
                numLinks = paths.getLong("num_paths");
                linkBytes = paths.getLong("path_bytes");

                JSONArray colors = graph.getJSONArray("colours");
                for (int i = 0; i < colors.length(); i++) {
                    JSONObject jColor = colors.getJSONObject(i);

                    CortexColor color = new CortexColor();
                    color.setSampleName(jColor.getString("sample"));
                    color.setTotalSequence(jColor.getLong("total_sequence"));
                    color.setTipClippingApplied(jColor.getBoolean("cleaned_tips"));
                    color.setLowCovgSupernodesRemoved(true);

                    this.colors.add(color);
                }
            }

            while ((line = buffered.readLine()) != null) {
                if (!line.isEmpty() && !line.startsWith("#")) {
                    break;
                }

                recordsStart += line.length() + 1;
                buffered.mark(100);
            }
        } catch (IOException e) {
            throw new IndianaException("Unable to load Cortex links file '" + linksFile.getAbsolutePath() + "'", e);
        }
    }

    private void moveToBeginningOfRecordsSection() {
        try {
            buffered.reset();
            recordsSeen = 0;
        } catch (IOException e) {
            try {
                buffered.close();

                InputStream fileStream = new FileInputStream(this.linksFile);
                InputStream gzipStream = new GZIPInputStream(fileStream);
                Reader decoder = new InputStreamReader(gzipStream, Charset.forName("UTF-8"));
                buffered = new BufferedReader(decoder);

                buffered.skip(recordsStart);
                recordsSeen = 0;
            } catch (IOException e2) {
                throw new IndianaException("Unable to restart iteration over CortexLinks file", e2);
            }
        }

        nextRecord = getNextRecord();
    }

    private CortexLinksRecord getNextRecord() {
        if (recordsSeen < numKmersWithLinks) {
            try {
                String[] kmerLine = buffered.readLine().split("\\s+");

                String kmer = kmerLine[0];
                int numLinks = Integer.valueOf(kmerLine[1]);

                List<CortexJunctionsRecord> cjs = new ArrayList<CortexJunctionsRecord>();

                for (int i = 0; i < numLinks; i++) {
                    String[] linkLine = buffered.readLine().split("\\s+");

                    String orientation = linkLine[0];
                    int numKmers = Integer.valueOf(linkLine[1]);
                    int numJunctions = Integer.valueOf(linkLine[2]);
                    int[] coverages = new int[numColors];

                    for (int c = 0; c < numColors; c++) {
                        coverages[c] = Integer.valueOf(linkLine[3 + c]);
                    }

                    String junctions = linkLine[linkLine.length - 1];

                    CortexJunctionsRecord cj = new CortexJunctionsRecord(orientation.equals("F"), numKmers, numJunctions, coverages, junctions);
                    cjs.add(cj);
                }

                recordsSeen++;

                return new CortexLinksRecord(kmer, cjs);
            } catch (IOException e) {
                throw new IndianaException("Unable to parse CortexLinks record", e);
            }
        }

        return null;
    }

    @Override
    public Iterator<CortexLinksRecord> iterator() {
        moveToBeginningOfRecordsSection();

        return this;
    }

    @Override
    public boolean hasNext() {
        return nextRecord != null;
    }

    @Override
    public CortexLinksRecord next() {
        CortexLinksRecord currentRecord = nextRecord;

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
            throw new IndianaException("Unable to close CortexLinks file", e);
        }
    }

    public String toString() {
        StringBuilder out = new StringBuilder();

        out.append("file: ").append(linksFile.getAbsolutePath()).append("\n");
        out.append("----\n");
        out.append("binary version: ").append(version).append("\n");
        out.append("kmer size: ").append(kmerSize).append("\n");
        out.append("link bytes: ").append(linkBytes).append("\n");
        out.append("colors: ").append(numColors).append("\n");
        for (int i = 0; i < colors.size(); i++) {
            CortexColor color = colors.get(i);

            out.append("-- Color ").append(i).append(" --\n");
            out.append("  sample: '").append(color.getSampleName()).append("'\n");
            out.append("  total_sequence: ").append(color.getTotalSequence()).append("\n");
            out.append("  cleaned_tips: ").append(color.isTipClippingApplied()).append("\n");
            out.append("  cleaned_supernodes: ").append(color.isLowCovgSupernodesRemoved()).append("\n");
        }
        out.append("----\n");
        out.append("kmers in graph: ").append(numKmersInGraph).append("\n");
        out.append("kmers with links: ").append(numKmersWithLinks).append("\n");
        out.append("paths: ").append(numLinks).append("\n");
        out.append("----\n");

        return out.toString();
    }

    public int getVersion() { return version; }
    public int getKmerSize() { return kmerSize; }
    public int getNumColors() { return numColors; }
    public long getNumKmersInGraph() { return numKmersInGraph; }
    public long getNumKmersWithLinks() { return numKmersWithLinks; }
    public long getNumLinks() { return numLinks; }
    public long getLinkBytes() { return linkBytes; }
    public File getCortexLinksFile() { return linksFile; }

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
