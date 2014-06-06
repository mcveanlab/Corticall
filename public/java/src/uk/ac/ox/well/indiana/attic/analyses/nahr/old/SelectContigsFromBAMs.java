package uk.ac.ox.well.indiana.attic.analyses.nahr.old;

import net.sf.samtools.*;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.util.*;

public class SelectContigsFromBAMs extends Module {
    @Argument(fullName="bam", shortName="bams", doc="BAM file")
    public ArrayList<SAMFileReader> BAMS;

    @Argument(fullName="id", shortName="id", doc="Contig IDs")
    public ArrayList<String> IDS;

    @Output
    public File out;

    @Override
    public void execute() {
        Map<String, Set<String>> ids = new HashMap<String, Set<String>>();
        for (String id : IDS) {
            String[] pieces = id.split("\\.");
            String sampleName = pieces[0];
            String contigName = pieces[1];

            if (!ids.containsKey(sampleName)) {
                ids.put(sampleName, new HashSet<String>());
            }

            ids.get(sampleName).add(contigName);
        }

        List<SAMRecord> contigs = new ArrayList<SAMRecord>();

        SAMFileHeader sfh = new SAMFileHeader();
        sfh.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        for (SAMSequenceRecord ssr : BAMS.get(0).getFileHeader().getSequenceDictionary().getSequences()) {
            sfh.addSequence(ssr);
        }

        for (SAMFileReader BAM : BAMS) {
            for (SAMReadGroupRecord rg : BAM.getFileHeader().getReadGroups()) {
                sfh.addReadGroup(rg);
            }

            for (SAMRecord contig : BAM) {
                String sampleName = contig.getReadGroup().getSample();
                String contigName = contig.getReadName();

                if (ids.containsKey(sampleName) && ids.get(sampleName).contains(contigName)) {
                    contigs.add(contig);
                }
            }
        }

        SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
        sfwf.setCreateIndex(true);
        SAMFileWriter sfw = sfwf.makeSAMOrBAMWriter(sfh, false, out);

        for (SAMRecord contig : contigs) {
            sfw.addAlignment(contig);
        }

        sfw.close();
    }
}
