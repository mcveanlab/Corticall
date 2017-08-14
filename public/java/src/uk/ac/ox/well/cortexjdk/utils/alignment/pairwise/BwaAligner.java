package uk.ac.ox.well.cortexjdk.utils.alignment.pairwise;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.ProcessExecutor;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class BwaAligner implements ExternalAligner {
    private final String bwaPath = System.getProperty("user.home") + "/bin/bwa";

    public List<SAMRecord> align(List<ReferenceSequence> queries, File targets) {
        try {
            File tempQuery = File.createTempFile("query", ".fa");

            PrintStream qw = new PrintStream(tempQuery);
            for (ReferenceSequence query : queries) {
                qw.println(">" + query.getName());
                qw.println(new String(query.getBases()));
            }
            qw.close();

            String result = ProcessExecutor.executeAndReturnResult(String.format("%s mem %s %s", bwaPath, targets.getAbsolutePath(), tempQuery.getAbsolutePath()));

            tempQuery.delete();

            List<SAMRecord> recs = new ArrayList<>();

            FastaSequenceFile fa = new FastaSequenceFile(targets, true);
            SAMFileHeader sfh = new SAMFileHeader();
            sfh.setSequenceDictionary(fa.getSequenceDictionary());
            sfh.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            for (String samLine : result.split("\\n")) {
                if (!samLine.isEmpty() && !samLine.startsWith("@")) {
                    //System.out.println("samLine: '" + samLine + "'");

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

            File tempTargetIdx = File.createTempFile("target", ".fa.fai");

            ProcessExecutor.executeAndReturnResult(String.format("%s index %s", bwaPath, tempTarget.getAbsolutePath()));
            String result = ProcessExecutor.executeAndReturnResult(String.format("%s mem %s %s", bwaPath, tempTarget.getAbsolutePath(), tempQuery.getAbsolutePath()));

            List<SAMRecord> recs = new ArrayList<>();

            SAMFileHeader sfh = new SAMFileHeader();
            sfh.setSortOrder(SAMFileHeader.SortOrder.unsorted);
            SAMSequenceDictionary ssd = new SAMSequenceDictionary();
            ssd.addSequence(new SAMSequenceRecord("target", target.length()));
            sfh.setSequenceDictionary(ssd);

            for (String samLine : result.split("\n")) {
                if (!samLine.isEmpty() && !samLine.startsWith("@")) {
                    System.out.println("test: '" + samLine + "'");

                    recs.add(new SAMLineParser(sfh).parseLine(samLine));
                }
            }

            tempQuery.delete();
            tempTarget.delete();
            tempTargetIdx.delete();

            return recs;
        } catch (IOException e) {
            throw new CortexJDKException("IOException: " + e);
        }
    }

    public List<SAMRecord> align(String query, File targets) {
        try {
            File tempQuery = File.createTempFile("query", ".fa");

            PrintStream qw = new PrintStream(tempQuery);
            qw.println(">query");
            qw.println(query);
            qw.close();

            String result = ProcessExecutor.executeAndReturnResult(String.format("%s mem %s %s", bwaPath, targets.getAbsolutePath(), tempQuery.getAbsolutePath()));

            List<SAMRecord> recs = new ArrayList<>();

            FastaSequenceFile fa = new FastaSequenceFile(targets, true);
            SAMFileHeader sfh = new SAMFileHeader();
            sfh.setSequenceDictionary(fa.getSequenceDictionary());
            sfh.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            for (String samLine : result.split("\n")) {
                if (!samLine.isEmpty() && !samLine.startsWith("@")) {
                    recs.add(new SAMLineParser(sfh).parseLine(samLine));
                }
            }

            tempQuery.delete();

            return recs;
        } catch (IOException e) {
            throw new CortexJDKException("IOException: " + e);
        }
    }
}
