package uk.ac.ox.well.indiana.utils.alignment.pairwise;

import com.google.common.base.Joiner;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.ProcessExecutor;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public class LastzAligner {
    //private final String lastzPath = "/Users/kiran/opt/lastz-distrib-1.03.66/bin/lastz_D";
    //private final String lastzPath = "lastz";
    //private final String lastzPath = "/Users/kiran/opt/lastz-distrib-1.03.66/bin/lastz";
    private final String lastzPath = System.getProperty("user.home") + "/opt/lastz-distrib-1.03.66/bin/lastz";

    public List<SAMRecord> align(String query, File targets) {
        try {
            File tempQuery = File.createTempFile("query", ".fa");

            PrintStream qw = new PrintStream(tempQuery);
            qw.println(">query");
            qw.println(query);
            qw.close();

            String hsx = targets.getAbsolutePath().replaceAll(".fasta$", ".hsx");
            String result = ProcessExecutor.executeAndReturnResult(String.format("%s %s[multiple] %s --format=%s --queryhspbest=1", lastzPath, hsx, tempQuery.getAbsolutePath(), "sam-"));

            tempQuery.delete();

            List<SAMRecord> recs = new ArrayList<SAMRecord>();

            FastaSequenceFile fa = new FastaSequenceFile(targets, true);
            SAMFileHeader sfh = new SAMFileHeader();
            sfh.setSequenceDictionary(fa.getSequenceDictionary());
            sfh.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            for (String samLine : result.split("\n")) {
                //System.out.println(samLine);

                if (!samLine.isEmpty()) {
                    recs.add(new SAMLineParser(sfh).parseLine(samLine));
                }
            }

            return recs;
        } catch (IOException e) {
            throw new IndianaException("IOException: " + e);
        }
    }

    public void align(ReferenceSequence query, File targets, String format) {
        try {
            File tempQuery = File.createTempFile("query", ".fa");

            PrintStream qw = new PrintStream(tempQuery);
            qw.println(">" + query.getName());
            qw.println(new String(query.getBases()));
            qw.close();

            String hsx = targets.getAbsolutePath().replaceAll(".fasta$", ".hsx");
            String result = ProcessExecutor.executeAndReturnResult(String.format("%s %s[multiple] %s --format=%s --queryhspbest=1", lastzPath, hsx, tempQuery.getAbsolutePath(), format));

            tempQuery.delete();

            System.out.println(result);
        } catch (IOException e) {
            throw new IndianaException("IOException: " + e);
        }
    }

    public Map<String, String[]> alignAll(Set<ReferenceSequence> queries, File targets) {
        try {
            File tempQueries = File.createTempFile("queries", ".fa");

            PrintStream qw = new PrintStream(tempQueries);
            for (ReferenceSequence query : queries) {
                qw.println(">" + query.getName());
                qw.println(new String(query.getBases()));
            }
            qw.close();

            String hsx = targets.getAbsolutePath().replaceAll(".fasta$", ".hsx");
            String result = ProcessExecutor.executeAndReturnResult(String.format("%s %s[multiple] %s --format=%s --queryhspbest=1", lastzPath, hsx, tempQueries.getAbsolutePath(), "sam-"));

            System.out.println(result);

            Map<String, Set<String[]>> alignments = new HashMap<String, Set<String[]>>();
            for (String line : result.split("\n")) {
                String[] fields = line.split("\\s+");

                String contigName = fields[0];

                if (!alignments.containsKey(contigName)) {
                    alignments.put(contigName, new HashSet<String[]>());
                }

                alignments.get(contigName).add(fields);
            }

            Map<String, String[]> results = new HashMap<String, String[]>();

            for (String contigName : alignments.keySet()) {
                Set<String> cigars = new HashSet<String>();

                for (String[] fields : alignments.get(contigName)) {
                    cigars.add(fields[5]);
                }

                if (cigars.size() == 1) {
                    String[] fields = alignments.get(contigName).iterator().next();
                    if (alignments.get(contigName).size() > 1) {
                        fields[4] = "0";
                    }

                    results.put(contigName, fields);
                }
            }

            tempQueries.delete();

            return results;
        } catch (IOException e) {
            throw new IndianaException("IOException: " + e);
        }
    }
}
