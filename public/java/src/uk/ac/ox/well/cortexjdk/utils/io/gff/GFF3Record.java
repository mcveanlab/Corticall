package uk.ac.ox.well.cortexjdk.utils.io.gff;

import com.google.common.base.Joiner;
import htsjdk.samtools.util.Interval;

import java.util.TreeMap;

public class GFF3Record implements Comparable<GFF3Record> {
    public enum Strand { POSITIVE, NEGATIVE }

    private String seqid;
    private String source;
    private String type;
    private int start;
    private int end;
    private String score;
    private Strand strand;
    private String phase;
    private TreeMap<String, String> attributes = new TreeMap<>();

    public GFF3Record() {}

    public GFF3Record(String line) {
        String[] fields = line.split("\t+");

        seqid = fields[0];
        source = fields[1];
        type = fields[2];
        start = Integer.valueOf(fields[3]);
        end = Integer.valueOf(fields[4]);
        score = fields[5];
        strand = fields[6].equals("+") ? Strand.POSITIVE : Strand.NEGATIVE;
        phase = fields[7];

        for (String attribute : fields[8].split(";")) {
            String[] keyvalue = attribute.split("=");

            attributes.put(keyvalue[0], keyvalue[1]);
        }
    }

    //public static GFF3Record construct(String seqid, String source, String type, int start, int end, int score, Strand strand, String phase)

    public GFF3Record(GFF3Record record) {
        this.seqid = record.seqid;
        this.source = record.source;
        this.type = record.type;
        this.start = record.start;
        this.end = record.end;
        this.score = record.score;
        this.strand = record.strand;
        this.phase = record.phase;
        this.attributes = record.attributes;
    }

    public String getSeqid()  { return seqid;  }
    public String getSource() { return source; }
    public String getType()   { return type;   }
    public int    getStart()  { return start;  }
    public int    getEnd()    { return end;    }
    public String getScore()  { return score;  }
    public Strand getStrand() { return strand; }
    public String getPhase()  { return phase;  }
    public boolean hasAttribute(String key) { return attributes.containsKey(key); }
    public String getAttribute(String key) { return attributes.get(key); }
    public TreeMap<String, String> getAttributes() { return attributes; }
    public Interval getInterval() { return new Interval(seqid, start, end); }

    public void setSeqid(String seqid) { this.seqid = seqid; }
    public void setSource(String source) { this.source = source; }
    public void setType(String type) { this.type = type; }
    public void setStart(int start) { this.start = start; }
    public void setEnd(int end) { this.end = end; }
    public void setScore(String score) { this.score = score; }
    public void setStrand(Strand strand) { this.strand = strand; }
    public void setPhase(String phase) { this.phase = phase; }
    public void setAttribute(String key, String value) { this.attributes.put(key, value); }
    public void setAttributes(TreeMap<String, String> attributes) { this.attributes = attributes; }

    public String toString() {
        return seqid  + "\t" +
               source + "\t" +
               type   + "\t" +
               start  + "\t" +
               end    + "\t" +
               score  + "\t" +
               (strand == Strand.POSITIVE ? "+" : "-") + "\t" +
               phase  + "\t" +
               Joiner.on(";").withKeyValueSeparator("=").join(attributes);
    }

    public int hashCode() {
        return this.toString().hashCode();
    }

    @Override
    public int compareTo(GFF3Record o) {
        if (this.getSeqid().equalsIgnoreCase(o.getSeqid())) {
            if (this.getStart() == o.getStart()) {
                if (this.getEnd() == o.getEnd()) {
                    return 0;
                } else {
                    return this.getEnd() < o.getEnd() ? -1 : 1;
                }
            } else {
                return this.getStart() <= o.getStart() ? -1 : 1;
            }
        }

        return this.getSeqid().compareTo(o.getSeqid());
    }
}
