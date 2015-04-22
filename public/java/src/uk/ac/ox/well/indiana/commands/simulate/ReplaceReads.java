package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

public class ReplaceReads extends Module {
    @Argument(fullName="from", shortName="f", doc="Replace reads from here")
    public SAMFileReader FROM;

    @Argument(fullName="with", shortName="w", doc="Replace reads with these")
    public SAMFileReader WITH;

    @Output
    public File out;

    @Override
    public void execute() {
        Map<String, Set<SAMRecord>> reads = new LinkedHashMap<String, Set<SAMRecord>>();
        Map<String, Set<SAMRecord>> replacements = new LinkedHashMap<String, Set<SAMRecord>>();

        log.info("Loading reads...");
        int numRecords = 0;
        for (SAMRecord fr : FROM) {
            if (!reads.containsKey(fr.getReadName())) {
                reads.put(fr.getReadName(), new HashSet<SAMRecord>());
            }

            reads.get(fr.getReadName()).add(fr);

            numRecords++;
        }

        int numReplacements = 0;
        for (SAMRecord wr : WITH) {
            if (!replacements.containsKey(wr.getReadName())) {
                replacements.put(wr.getReadName(), new HashSet<SAMRecord>());
            }

            wr.setHeader(FROM.getFileHeader());
            wr.setReferenceIndex(FROM.getFileHeader().getSequenceIndex(wr.getReferenceName()));
            wr.setAttribute("RG", FROM.getFileHeader().getReadGroups().get(0).getReadGroupId());

            replacements.get(wr.getReadName()).add(wr);

            numReplacements++;
        }
        log.info("  {} original records", numRecords);
        log.info("  {} replacements", numReplacements);

        log.info("Writing reads...");
        //SAMFileWriter sfw = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(FROM.getFileHeader(), false, out);
        SAMFileWriter sfw = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(FROM.getFileHeader(), false, out);
        int numWritten = 0;

        for (String readName : reads.keySet()) {
            Set<SAMRecord> records = (replacements.containsKey(readName) ? replacements.get(readName) : reads.get(readName));

            for (SAMRecord read : records) {
                /*
                if (read.getReadGroup() == null) {
                    read.setAttribute("RG", FROM.getFileHeader().getReadGroups().get(0).getReadGroupId());
                }
                */

                sfw.addAlignment(read);

                numWritten++;
            }
        }

        sfw.close();

        log.info("  {} records written", numWritten);
    }
}
