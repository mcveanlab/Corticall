package uk.ac.ox.well.indiana.utils.alignment.pairwise;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMLineParser;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.ProcessExecutor;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class BwaAligner {
    private final String bwaPath = System.getProperty("user.home") + "/bin/bwa";

    public List<SAMRecord> align(String query, File targets) {
        try {
            File tempQuery = File.createTempFile("query", ".fa");

            PrintStream qw = new PrintStream(tempQuery);
            qw.println(">query");
            qw.println(query);
            qw.close();

            String result = ProcessExecutor.executeAndReturnResult(String.format("%s mem %s %s", bwaPath, targets.getAbsolutePath(), tempQuery.getAbsolutePath()));

            List<SAMRecord> recs = new ArrayList<SAMRecord>();

            FastaSequenceFile fa = new FastaSequenceFile(targets, true);
            SAMFileHeader sfh = new SAMFileHeader();
            sfh.setSequenceDictionary(fa.getSequenceDictionary());
            sfh.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            for (String samLine : result.split("\n")) {
                recs.add(new SAMLineParser(sfh).parseLine(samLine));
            }

            tempQuery.delete();

            return recs;
        } catch (IOException e) {
            throw new IndianaException("IOException: " + e);
        }
    }
}
