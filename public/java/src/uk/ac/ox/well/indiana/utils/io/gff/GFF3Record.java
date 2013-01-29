package uk.ac.ox.well.indiana.utils.io.gff;

import com.google.common.base.Joiner;

import java.util.TreeMap;

public class GFF3Record {
    public enum Strand { POSITIVE, NEGATIVE }

    private String seqid;
    private String source;
    private String type;
    private int start;
    private int end;
    private String score;
    private Strand strand;
    private String phase;
    private TreeMap<String, String> attributes = new TreeMap<String, String>();

    public GFF3Record(String line) {
        String[] fields = line.split("\t");

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

    public String getSeqid()  { return seqid;  }
    public String getSource() { return source; }
    public String getType()   { return type;   }
    public int    getStart()  { return start;  }
    public int    getEnd()    { return end;    }
    public String getScore()  { return score;  }
    public Strand getStrand() { return strand; }
    public String getPhase()  { return phase;  }
    public TreeMap<String, String> getAttributes() { return attributes; }

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
}
