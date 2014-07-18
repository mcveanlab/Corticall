package uk.ac.ox.well.indiana.attic.analyses.nahr;

import htsjdk.samtools.*;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class MergeContigsAndOverlappingReads extends Module {
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

        SAMFileHeader sfh = CONTIGS.getFileHeader();
        for (SAMReadGroupRecord srg : READS.getFileHeader().getReadGroups()) {
            sfh.addReadGroup(srg);
        }

        SAMFileWriter sfw = sfwf.makeBAMWriter(sfh, false, out);

        for (SAMRecord contig : CONTIGS) {
            if (CONTIG_NAMES.contains(contig.getReadName())) {
                contig.setAttribute("XC", 0);
                sfw.addAlignment(contig);

                int alignmentStart = contig.getAlignmentStart();
                int alignmentEnd = contig.getAlignmentEnd();

                CigarElement firstCe = contig.getCigar().getCigarElements().get(0);
                CigarElement lastCe = contig.getCigar().getCigarElements().get(contig.getCigar().getCigarElements().size() - 1);

                if (firstCe.getOperator().equals(CigarOperator.S)) {
                    alignmentStart = contig.getAlignmentStart() - firstCe.getLength();
                }

                if (lastCe.getOperator().equals(CigarOperator.S)) {
                    alignmentEnd = contig.getAlignmentEnd() + lastCe.getLength();
                }

                //SAMRecordIterator sri = READS.queryOverlapping(contig.getReferenceName(), contig.getAlignmentStart(), contig.getAlignmentEnd());
                SAMRecordIterator sri = READS.queryOverlapping(contig.getReferenceName(), alignmentStart, alignmentEnd);
                while (sri.hasNext()) {
                    SAMRecord read = sri.next();

                    contig.setAttribute("XC", 1);

                    sfw.addAlignment(read);
                }
                sri.close();
            }
        }

        sfw.close();
    }
}
