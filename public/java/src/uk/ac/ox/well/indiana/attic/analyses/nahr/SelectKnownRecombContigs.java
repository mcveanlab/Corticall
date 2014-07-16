package uk.ac.ox.well.indiana.attic.analyses.nahr;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.File;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

public class SelectKnownRecombContigs extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public SAMFileReader BAM;

    @Argument(fullName="recombTable", shortName="rt", doc="Recomb event table")
    public File RECOMB_TABLE;

    @Output
    public File out;

    @Output(fullName="stats_out", shortName="so", doc="Stats output file")
    public PrintStream sout;

    @Override
    public void execute() {
        Set<String> partialSampleNames = new HashSet<String>();
        for (SAMReadGroupRecord rgr : BAM.getFileHeader().getReadGroups()) {
            String sampleName = rgr.getSample();
            String readGroupId = rgr.getReadGroupId();
            readGroupId = readGroupId.replaceAll("\\..+", "");

            String partialSampleName = sampleName + "/" + readGroupId;
            partialSampleNames.add(partialSampleName);
        }

        TableReader tr = new TableReader(RECOMB_TABLE);
        Set<Map<String, String>> recombEvents = new HashSet<Map<String, String>>();

        for (Map<String, String> te : tr) {
            String fullSampleName = te.get("sample");
            for (String partialSampleName : partialSampleNames) {
                if (fullSampleName.contains(partialSampleName)) {
                    recombEvents.add(te);
                }
            }
        }

        Set<SAMRecord> spanningContigs = new HashSet<SAMRecord>();

        TableWriter tw = new TableWriter(sout);

        for (Map<String, String> te : recombEvents) {
            String chrom = te.get("chrom");
            int start = Integer.valueOf(te.get("co_pos_min"));
            int end = Integer.valueOf(te.get("co_pos_max"));

            Interval fullInterval = new Interval(chrom, start, end);

            int contigsFound = 0;
            SAMRecordIterator sri = BAM.queryOverlapping(chrom, start, end);
            while (sri.hasNext()) {
                SAMRecord contig = sri.next();

                Interval contigInterval = new Interval(contig.getReferenceName(), contig.getAlignmentStart(), contig.getAlignmentEnd());

                if (fullInterval.getIntersectionLength(contigInterval) == fullInterval.length()) {
                    spanningContigs.add(contig);

                    log.info("chrom={} start={} end={} contig={}", chrom, start, end, contig.getSAMString());

                    contigsFound++;
                }

            }
            sri.close();

            Map<String, String> entry = new LinkedHashMap<String, String>();
            entry.put("sample", te.get("sample"));
            entry.put("chrom", te.get("chrom"));
            entry.put("co_pos_min", te.get("co_pos_min"));
            entry.put("co_pos_max", te.get("co_pos_max"));
            entry.put("contigs_found", String.valueOf(contigsFound));

            tw.addEntry(entry);
        }

        SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
        sfwf.setCreateIndex(true);

        SAMFileWriter sfw = sfwf.makeBAMWriter(BAM.getFileHeader(), false, out);

        for (SAMRecord spanningContig : spanningContigs) {
            sfw.addAlignment(spanningContig);
        }

        sfw.close();
    }
}
