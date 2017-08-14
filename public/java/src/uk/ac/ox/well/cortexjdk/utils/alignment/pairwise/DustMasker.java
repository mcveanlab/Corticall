package uk.ac.ox.well.cortexjdk.utils.alignment.pairwise;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.ProcessExecutor;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public class DustMasker {
    private File dustmaskerPath = new File("/Users/kiran/opt/ncbi-blast-2.5.0+/bin/dustmasker");

    public DustMasker() {
        if (!dustmaskerPath.exists()) {
            dustmaskerPath = new File(ProcessExecutor.executeAndReturnResult("which dustmasker").trim());

            if (!dustmaskerPath.exists()) {
                throw new CortexJDKException("'dustmasker' is not on the PATH.");
            }
        }
    }

    public Map<Integer, Map<String, Object>> mask(Collection<ReferenceSequence> contigs) {
        try {
            Map<Integer, Map<String, Object>> ms = new HashMap<>();

            File tempQueries = File.createTempFile("query", ".fa");

            PrintStream qw = new PrintStream(tempQueries);
            for (ReferenceSequence rseq : contigs) {
                qw.println(">" + rseq.getName());
                qw.println(rseq.getBaseString());

                Map<String, Object> m = new HashMap<>();
                m.put("dust", false);

                ms.put(Integer.valueOf(rseq.getName()), m);
            }
            qw.close();

            String cmd = String.format(
                    "%s -in %s -outfmt acclist",
                    dustmaskerPath.getAbsolutePath(),
                    tempQueries.getAbsolutePath()
            );

            String result = ProcessExecutor.executeAndReturnResult(cmd);

            for (String line : result.split("\n")) {
                String[] pieces = line.split("\\s+");

                int contigNum = Integer.valueOf(pieces[0].replaceAll("^>", ""));
                ms.get(contigNum).put("dust", true);
            }

            return ms;
        } catch (IOException e) {
            return null;
        }
    }
}
