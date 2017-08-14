package uk.ac.ox.well.cortexjdk.utils.alignment.pairwise;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.ProcessExecutor;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public class LastzAligner implements ExternalAligner {
    private final String lastzDir = System.getProperty("user.home") + "/opt/lastz-distrib-1.03.66";
    private final String lastzPath = lastzDir + "/bin/lastz";
    private final String hsxPath = lastzDir + "/tools/build_fasta_hsx.py";

    public List<SAMRecord> align(String query, File targets) {
        try {
            File tempQueries = File.createTempFile("query", ".fa");

            PrintStream qw = new PrintStream(tempQueries);
            qw.println(">query");
            qw.println(query);
            qw.close();

            String hsx = targets.getAbsolutePath().replaceAll(".fasta$", ".hsx");
            String result = ProcessExecutor.executeAndReturnResult(String.format("%s %s[multiple] %s --format=%s --queryhspbest=1", lastzPath, hsx, tempQueries.getAbsolutePath(), "sam-"));

            tempQueries.delete();

            List<SAMRecord> recs = new ArrayList<>();

            FastaSequenceFile fa = new FastaSequenceFile(targets, true);
            SAMFileHeader sfh = new SAMFileHeader();
            sfh.setSequenceDictionary(fa.getSequenceDictionary());
            sfh.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            for (String samLine : result.split("\n")) {
                if (!samLine.isEmpty()) {
                    recs.add(new SAMLineParser(sfh).parseLine(samLine));
                }
            }

            return recs;
        } catch (IOException e) {
            throw new CortexJDKException("IOException: " + e);
        }
    }

    public List<SAMRecord> align(List<ReferenceSequence> queries, File targets) {
        try {
            File tempQuery = File.createTempFile("query", ".fa");

            PrintStream qw = new PrintStream(tempQuery);
            for (ReferenceSequence query : queries) {
                qw.println(">" + query.getName());
                qw.println(new String(query.getBases()));
            }
            qw.close();

            String hsx = targets.getAbsolutePath().replaceAll(".fasta$", ".hsx");
            if (!new File(hsx).exists()) {
                ProcessExecutor.execute(String.format("%s %s > %s", hsxPath, targets.getAbsolutePath(), hsx));
            }

            String result = ProcessExecutor.executeAndReturnResult(String.format("%s %s[multiple] %s --format=%s --queryhspbest=1 --ambiguous=iupac", lastzPath, hsx, tempQuery.getAbsolutePath(), "sam-"));

            tempQuery.delete();

            List<SAMRecord> recs = new ArrayList<>();

            FastaSequenceFile fa = new FastaSequenceFile(targets, true);
            SAMFileHeader sfh = new SAMFileHeader();
            sfh.setSequenceDictionary(fa.getSequenceDictionary());
            sfh.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            for (String samLine : result.split("\n")) {
                if (!samLine.isEmpty()) {
                    recs.add(new SAMLineParser(sfh).parseLine(samLine));
                }
            }

            return recs;
        } catch (IOException e) {
            throw new CortexJDKException("IOException: " + e);
        }
    }

    public List<SAMRecord> align(String query, String target) {
        try {
            File tempQuery = File.createTempFile("query", ".fa");

            PrintStream qw = new PrintStream(tempQuery);
            qw.println(">query");
            qw.println(query);
            qw.close();

            File tempTarget = File.createTempFile("target", ".fa");

            PrintStream tw = new PrintStream(tempTarget);
            tw.println(">target");
            tw.println(target);
            tw.close();

            //String result = ProcessExecutor.executeAndReturnResult(String.format("%s mem %s %s", bwaPath, tempTarget.getAbsolutePath(), tempQuery.getAbsolutePath()));
            //String hsx = targets.getAbsolutePath().replaceAll(".fasta$", ".hsx");
            //String result = ProcessExecutor.executeAndReturnResult(String.format("%s %s[multiple] %s --format=%s --queryhspbest=1", lastzPath, tempTarget.getAbsolutePath(), tempQuery.getAbsolutePath(), "sam-"));
            String result = ProcessExecutor.executeAndReturnResult(String.format("%s %s %s --notransition --step=20 --nogapped --format=%s --queryhspbest=1", lastzPath, tempTarget.getAbsolutePath(), tempQuery.getAbsolutePath(), "sam-"));

            List<SAMRecord> recs = new ArrayList<>();

            SAMFileHeader sfh = new SAMFileHeader();
            sfh.setSortOrder(SAMFileHeader.SortOrder.unsorted);
            SAMSequenceDictionary ssd = new SAMSequenceDictionary();
            ssd.addSequence(new SAMSequenceRecord("target", target.length()));
            sfh.setSequenceDictionary(ssd);

            for (String samLine : result.split("\n")) {
                if (!samLine.isEmpty() && !samLine.startsWith("@")) {
                    recs.add(new SAMLineParser(sfh).parseLine(samLine));
                }
            }

            tempQuery.delete();
            tempTarget.delete();

            return recs;
        } catch (IOException e) {
            throw new CortexJDKException("IOException: " + e);
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
            throw new CortexJDKException("IOException: " + e);
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

            Map<String, Set<String[]>> alignments = new HashMap<>();
            for (String line : result.split("\n")) {
                String[] fields = line.split("\\s+");

                String contigName = fields[0];

                if (!alignments.containsKey(contigName)) {
                    alignments.put(contigName, new HashSet<>());
                }

                alignments.get(contigName).add(fields);
            }

            Map<String, String[]> results = new HashMap<>();

            for (String contigName : alignments.keySet()) {
                Set<String> cigars = new HashSet<>();

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
            throw new CortexJDKException("IOException: " + e);
        }
    }
}
