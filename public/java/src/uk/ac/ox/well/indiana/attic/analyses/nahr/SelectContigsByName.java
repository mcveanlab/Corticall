package uk.ac.ox.well.indiana.attic.analyses.nahr;

import net.sf.samtools.*;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.util.HashSet;

public class SelectContigsByName extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (BAM)")
    public SAMFileReader CONTIGS;

    @Argument(fullName="contigNames", shortName="cn", doc="Contig names")
    public HashSet<String> CONTIG_NAMES;

    @Argument(fullName="reads", shortName="r", doc="Reads (BAM)")
    public SAMFileReader READS;

    @Output
    public File out;

    @Override
    public void execute() {
        SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
        sfwf.setCreateIndex(true);
        SAMFileWriter sfw = sfwf.makeBAMWriter(CONTIGS.getFileHeader(), false, out);

        SAMFileHeader sfh = CONTIGS.getFileHeader();
        for (SAMReadGroupRecord srg : READS.getFileHeader().getReadGroups()) {
            sfh.addReadGroup(srg);
        }

        for (SAMRecord contig : CONTIGS) {
            if (CONTIG_NAMES.contains(contig.getReadName())) {
                sfw.addAlignment(contig);

                SAMRecord read;
                while ((read = READS.queryOverlapping(contig.getReferenceName(), contig.getAlignmentStart(), contig.getAlignmentEnd()).next()) != null) {
                    sfw.addAlignment(read);
                }
            }
        }

        sfw.close();
    }
}
