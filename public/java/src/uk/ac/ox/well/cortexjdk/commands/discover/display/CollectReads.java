package uk.ac.ox.well.cortexjdk.commands.discover.display;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;

import java.io.PrintStream;
import java.util.*;

public class CollectReads extends Module {
    @Argument(fullName="contigReads", shortName="cr", doc="Contig reads")
    public SAMFileReader CONTIG_READS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, SAMRecord> readsEnd1 = new HashMap<>();
        Map<String, SAMRecord> readsEnd2 = new HashMap<>();

        log.info("Loading contig reads:");
        for (SAMRecord cr : CONTIG_READS) {
            if (cr.getMappingQuality() > 0) {
                Map<String, SAMRecord> reads = cr.getFirstOfPairFlag() ? readsEnd1 : readsEnd2;

                reads.put(cr.getReadName(), cr);
            }
        }

        log.info("Printing table:");
        for (String sn : readsEnd1.keySet()) {
            if (readsEnd1.containsKey(sn) && readsEnd2.containsKey(sn)) {
                SAMRecord sr1 = readsEnd1.get(sn);
                SAMRecord sr2 = readsEnd2.get(sn);

                out.println(Joiner.on("\t").join(
                        sr1.getReadName(), sr1.getContig(), sr1.getStart(), sr1.getEnd(), sr1.getMappingQuality(),
                                           sr2.getContig(), sr2.getStart(), sr2.getEnd(), sr2.getMappingQuality()
                ));
            }
        }
    }
}
