package uk.ac.ox.well.cortexjdk.commands.discover.display;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;

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

        log.info("Loading contig reads:");
        for (SAMRecord cr : CONTIG_READS) {
            if (cr.getMappingQuality() > 0) {
                Map<String, SAMRecord> reads = cr.getFirstOfPairFlag() ? readsEnd1 : readsEnd2;
                reads.put(cr.getReadName(), cr);
            }
        }
        log.info("  loaded {} reads [{}, {}]", readsEnd1.size() + readsEnd2.size(), readsEnd1.size(), readsEnd2.size());

        for (SAMFileReader sfr : Arrays.asList(REF0_READS, REF1_READS)) {
            ProgressMeter pm = new ProgressMeterFactory()
                    .header("Loading ref reads...")
                    .message("ref reads processed")
                    .updateRecord(1000000)
                    .make(log);

            for (SAMRecord cr : sfr) {
                Map<String, SAMRecord> reads = cr.getFirstOfPairFlag() ? readsEnd1 : readsEnd2;
                Map<String, Set<String>> chrs = cr.getFirstOfPairFlag() ? chrsEnd1 : chrsEnd2;

                if (cr.getMappingQuality() > 0 && reads.containsKey(cr.getReadName()) && cr.getContig() != null) {
                    if (!chrs.containsKey(cr.getReadName())) {
                        chrs.put(cr.getReadName(), new HashSet<>());
                    }

                    chrs.get(cr.getReadName()).add(cr.getContig());
                }

                pm.update();
            }
        }

        for (String sn : readsEnd1.keySet()) {
            if (readsEnd1.containsKey(sn) && readsEnd2.containsKey(sn) && chrsEnd1.containsKey(sn) && chrsEnd2.containsKey(sn)) {
                SAMRecord cr1 = readsEnd1.get(sn);
                SAMRecord cr2 = readsEnd2.get(sn);

                String chrs1 = Joiner.on(",").join(chrsEnd1.get(sn));
                String chrs2 = Joiner.on(",").join(chrsEnd2.get(sn));

                out.println(Joiner.on("\t").join(
                    cr1.getReadName(), cr1.getContig(), cr1.getStart(), cr1.getEnd(), cr1.getMappingQuality(), cr1.getIntegerAttribute("NM"), chrs1,
                    cr2.getReadName(), cr2.getContig(), cr2.getStart(), cr2.getEnd(), cr2.getMappingQuality(), cr2.getIntegerAttribute("NM"), chrs2
                ));
            }
        }
    }
}
