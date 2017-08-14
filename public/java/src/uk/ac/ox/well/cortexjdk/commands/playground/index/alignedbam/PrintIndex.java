package uk.ac.ox.well.cortexjdk.commands.playground.index.alignedbam;

import com.google.common.base.Joiner;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Interval;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexBinaryKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;

import java.io.File;
import java.io.PrintStream;
import java.util.Set;
import java.util.TreeSet;

public class PrintIndex extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM")
    public File SAM_FILE;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        SamReader sr = SamReaderFactory.make()
                .setOption(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, true)
                .setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, true)
                .open(SAM_FILE);

        KmerIndex ki = new KmerIndex(SAM_FILE, KMER_SIZE, false);

        log.info("recs: {}", ki.getNumRecords());

        for (int i = 0; i < ki.getNumRecords(); i++) {
            Pair<CortexBinaryKmer, long[]> p = ki.getRecord(i);

            SAMFileSpan sfs = new BAMFileSpan(new Chunk(p.getSecond()[0], p.getSecond()[1]));

            Set<String> intervals = new TreeSet<>();

            SAMRecordIterator recs = sr.indexing().iterator(sfs);

            while (recs.hasNext()) {
                SAMRecord srecord = recs.next();
                Interval interval = new Interval(srecord.getReferenceName(), srecord.getAlignmentStart(), srecord.getAlignmentEnd());

                String sinterval = String.format("%s:%d-%d", interval.getContig(), interval.getStart(), interval.getEnd());

                intervals.add(sinterval);
            }

            recs.close();

            log.info("{}: {} {}",
                    i,
                    new String(CortexRecord.decodeBinaryKmer(p.getFirst().getBinaryKmer(), KMER_SIZE, CortexRecord.getKmerBits(KMER_SIZE))),
                    Joiner.on(";").join(intervals)
            );
        }
    }
}
