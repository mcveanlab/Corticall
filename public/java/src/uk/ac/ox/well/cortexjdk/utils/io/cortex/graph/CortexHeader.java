package uk.ac.ox.well.cortexjdk.utils.io.cortex.graph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CortexHeader {
    private int version;
    private int kmerSize;
    private int kmerBits;
    private int numColors;

    private List<CortexColor> colors = new ArrayList<>();
    private Map<String, Integer> colorForSampleNamesMap = new HashMap<>();
    private Map<Integer, String> sampleNameForColorMap = new HashMap<>();

    public int getColorForSampleName(String sampleName) { return colorForSampleNamesMap.get(sampleName); }
    public String getSampleNameForColor(int color) { return sampleNameForColorMap.get(color); }

    public int getVersion() { return version; }
    public int getKmerSize() { return kmerSize; }
    public int getKmerBits() { return kmerBits; }
    public int getNumColors() { return numColors; }
    public List<CortexColor> getColors() { return colors; }
    public CortexColor getColor(int c) { return colors.get(c); }
    public boolean hasColor(int color) { return (color < colors.size()); }

    public void setVersion(int version) { this.version = version; }
    public void setKmerSize(int kmerSize) { this.kmerSize = kmerSize; }
    public void setKmerBits(int kmerBits) { this.kmerBits = kmerBits; }
    public void setNumColors(int numColors) { this.numColors = numColors; }
    public void addColor(CortexColor color) {
        colors.add(color);

        colorForSampleNamesMap.put(color.getSampleName(), colors.size() - 1);
        sampleNameForColorMap.put(colors.size() - 1, color.getSampleName());
    }
}
