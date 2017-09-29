package uk.ac.ox.well.cortexjdk.utils.alignment.pairwise;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.ProcessExecutor;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

public class BlastAligner {
    private File blastPath = new File("/Users/kiran/opt/ncbi-blast-2.5.0+/bin/blastn");
    private File blastDbPath = new File("/Users/kiran/opt/ncbi-blast-2.5.0+/db");

    public BlastAligner() {
        if (!blastPath.exists()) {
            blastPath = new File(ProcessExecutor.executeAndReturnResult("which blastn").trim());

            if (!blastPath.exists()) {
                throw new CortexJDKException("'blastn' is not on the PATH.");
            }
        }

        if (!blastDbPath.exists()) {
            blastDbPath = blastPath.toPath().getParent().getParent().resolve("db").toFile();

            if (!blastDbPath.exists()) {
                throw new CortexJDKException("BLAST database directory is not in expected location.");
            }
        }
    }

    public Map<Integer, Map<String, Object>> align(Collection<ReferenceSequence> contigs) {
        try {
            File tempQueries = File.createTempFile("query", ".fa");

            PrintStream qw = new PrintStream(tempQueries);
            for (ReferenceSequence rseq : contigs) {
                qw.println(">" + rseq.getName());
                qw.println(rseq.getBaseString());
            }
            qw.close();

            String cmd = String.format(
                    "%s -db %s/vector -query %s -num_alignments 1 -outfmt '%s'",
                    blastPath.getAbsolutePath(),
                    blastDbPath.getAbsolutePath(),
                    tempQueries.getAbsolutePath(),
                    "6 qseqid sseqid qlen slen qseq sseq mismatch gaps evalue score stitle"
            );

            File tempCmd = File.createTempFile("cmd", ".sh");
            tempCmd.deleteOnExit();

            PrintStream cw = new PrintStream(tempCmd);
            cw.println(cmd);

            String result = ProcessExecutor.executeAndReturnResult("sh " + tempCmd);

            Map<Integer, Map<String, Object>> ls = new HashMap<>();

            for (String line : result.split("\n")) {
                String[] pieces = line.split("\\s+");

                if (pieces.length == 11) {
                    Map<String, Object> m = new HashMap<>();
                    m.put("qseqid", pieces[0]);
                    m.put("sseqid", pieces[1]);
                    m.put("qlen", pieces[2]);
                    m.put("slen", pieces[3]);
                    m.put("qseq", pieces[4]);
                    m.put("sseq", pieces[5]);
                    m.put("mismatch", pieces[6]);
                    m.put("gaps", pieces[7]);
                    m.put("evalue", pieces[8]);
                    m.put("score", pieces[9]);
                    m.put("stitle", pieces[10]);

                    //ls.add(m);

                    ls.put(Integer.valueOf(pieces[0]), m);
                }
            }

            return ls;
        } catch (IOException e) {
            return null;
        }
    }
}
