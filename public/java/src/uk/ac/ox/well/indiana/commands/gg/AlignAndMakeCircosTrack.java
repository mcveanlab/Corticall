package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.BwaAligner;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.ExternalAligner;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class AlignAndMakeCircosTrack extends Module {
    @Argument(fullName="query", shortName="q", doc="Query assembly")
    public FastaSequenceFile QUERY;

    @Argument(fullName="target", shortName="t", doc="Target assembly")
    public File TARGET;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        List<ReferenceSequence> queries = new ArrayList<ReferenceSequence>();

        ReferenceSequence rseq;
        while ((rseq = QUERY.nextSequence()) != null) {
            queries.add(rseq);
        }

        log.info("Aligning {} contigs...", queries.size());
        ExternalAligner ea = new BwaAligner();
        List<SAMRecord> alignments = ea.align(queries, TARGET);

        Map<String, Integer> alignmentCount = new HashMap<String, Integer>();
        for (SAMRecord alignment : alignments) {
            ContainerUtils.increment(alignmentCount, alignment.getReadName());
        }

        log.info("Writing alignment track...");
        int unalignedPos = 1;
        int mapped = 0, unmapped = 0, multiple = 0;

        for (SAMRecord alignment : alignments) {
            if (alignmentCount.get(alignment.getReadName()) == 1) {
                if (alignment.getReadUnmappedFlag()) {
                    out.printf("NA %d %d # %s\n", unalignedPos, unalignedPos + alignment.getReadLength(), alignment.getReadName());
                    unalignedPos += alignment.getReadLength();

                    unmapped++;
                } else {
                    out.printf("%s %d %d # %s\n", alignment.getReferenceName(), alignment.getAlignmentStart(), alignment.getAlignmentEnd(), alignment.getReadName());

                    mapped++;
                }
            } else {
                multiple++;
            }
        }

        log.info("  mapped={} unmapped={} multiple={} total={}", mapped, unmapped, multiple, queries.size());
    }
}
