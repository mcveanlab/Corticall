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

    @Argument(fullName="ref0Reads", shortName="r0", doc="Ref 0 reads")
    public SAMFileReader REF0_READS;

    @Argument(fullName="ref1Reads", shortName="r1", doc="Ref 1 reads")
    public SAMFileReader REF1_READS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, SAMRecord> readsEnd1 = new HashMap<>();
        Map<String, SAMRecord> readsEnd2 = new HashMap<>();
        Map<String, Set<String>> chrsEnd1 = new HashMap<>();
        Map<String, Set<String>> chrsEnd2 = new HashMap<>();

        for (SAMRecord cr : CONTIG_READS) {
            if (cr.getMappingQuality() > 0) {
                Map<String, SAMRecord> reads = cr.getFirstOfPairFlag() ? readsEnd1 : readsEnd2;
                reads.put(cr.getReadName(), cr);
            }
        }

        for (SAMFileReader sfr : Arrays.asList(REF0_READS, REF1_READS)) {
            for (SAMRecord cr : sfr) {
                Map<String, SAMRecord> reads = cr.getFirstOfPairFlag() ? readsEnd1 : readsEnd2;
                Map<String, Set<String>> chrs = cr.getFirstOfPairFlag() ? chrsEnd1 : chrsEnd2;

                if (cr.getMappingQuality() > 0 && reads.containsKey(cr.getReadName())) {
                    if (!chrs.containsKey(cr.getReadName())) {
                        chrs.put(cr.getReadName(), new HashSet<>());
                    }

                    chrs.get(cr.getReadName()).add(cr.getContig());
                }
            }
        }

        for (String sn : readsEnd1.keySet()) {
            if (readsEnd1.containsKey(sn) && readsEnd2.containsKey(sn)) {
                SAMRecord cr1 = readsEnd1.get(sn);
                SAMRecord cr2 = readsEnd2.get(sn);

                String chrs1 = Joiner.on(",").join(chrsEnd1.get(sn));
                String chrs2 = Joiner.on(",").join(chrsEnd2.get(sn));

                out.println(Joiner.on("\t").join(
                    cr1.getReadName(), cr1.getContig(), cr1.getStart(), cr1.getEnd(), chrs1,
                    cr2.getReadName(), cr2.getContig(), cr2.getStart(), cr2.getEnd(), chrs2
                ));
            }
        }
    }
}
